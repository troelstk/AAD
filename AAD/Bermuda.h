//
//  Bermuda.h
//  AAD
//
//  Created by Troels Tang Karlsen on 08/11/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef Bermuda_h
#define Bermuda_h

#include <cmath>
#include "utilities.h"
#include "gaussians.h"
#include "mrg32k3a.h"
#include "/usr/local/Cellar/armadillo/9.700.2/include/armadillo"

using namespace std;


template<class T> T LMM_BermudaSwaption(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF, vector<double> & exTimes,
                             double t, double Ta, double Tb, T r_fix, double notional,
                             int seed1, int seed2, int nPaths, int nPaths_presim, int nSteps_y, double yearly_payments, int dim_n)
{
    clock_t begin_time = clock();
    T tau(1.0/yearly_payments );
    int int_Ta = (int)(Ta), int_Tb = (int)(Tb);
    
    int no_basis_funcs = 3;
    mat beta(no_basis_funcs, exTimes.size() - 1);
    double lambda = 0.0; // Tikhonov parameter
    
    size_t M1 = initF.size();
    
    vec means(M1, fill::zeros);
    //mat C(M1, M1, fill::zeros);
    mat C(M1, M1, fill::zeros);
    
    // Fill covariance matrix
    for(int i = 0; i<M1; ++i){
        for(int j = 0; j<M1; ++j){
            C(i,j) = corr[i][j];
        }
    }
    mat corr2 = C( span(int_Ta,M1-1), span(int_Ta,M1-1) );
    // Find eigenvalues and vectors of correlation matrix
    mat P;
    vec eigval;
    eig_sym(eigval, P, corr2);
    // H has diagonal with eigenvalues
    mat H = diagmat(eigval);
    // Check if eigval is sorted ascending
    if( eigval.is_sorted() == 0 ) return 0.0;
    vec largest_n_eigvals = eigval.tail_rows(dim_n);
    // H2 is H but only largest n eigenvalues
    mat H_n = diagmat(largest_n_eigvals);
    // Lambda_n = sqrt(H_n)
    mat L_n = sqrt(H_n);
    // P_n is columns corresponding to largest n eigenvalue
    mat P_n = P.tail_cols(dim_n);
    
    // B_n is used to simulate from
    mat B_n = P_n * L_n;
    
    arma_rng::set_seed(seed1);
    
    vector<T> F;
    //mrg32k3a myRNG(seed1, seed2) ;
    //myRNG.init(dim_n);
    print("Init took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    
    
    
    //vector<double> gauss(dim_n);
    { // Scope to destruct everything declared in pre-simulation and first backward loop
        vector<vector<double>> SR(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> Libor(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> swap_vals(nPaths_presim, vector<double>(exTimes.size()));
        
        // Simulate swap rates, pre-simulation
        for(int t = 1; t<exTimes.size(); ++t ){
            //print("exTime is ", exTimes[t], " t is ", t);
            T swap_val(0.0);
            int int_ExTime = (int)(exTimes[t]);
            int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
            T dt( 1.0 /double(nSteps_y));
            T sqDt = sqrt(dt);

            for(int i = 0; i<nPaths_presim; ++i){
                mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
                vector<T> lnF = log(initF);
                for(int n = 0; n<nSteps; ++n) {
                    //print_DEBUG("Time is ", exTimes[t-1] + dt*(n+1));
                    for(int k = int_ExTime; k < int_Tb; ++k)
                    {
                        T sum(0.0);
                        for(int j = int_ExTime; j <= k; ++j ) {
                            double Fj = exp(lnF[j]);
                            sum += corr[k][j] * vol[j][j] * Fj / (1.0 + tau * Fj);
                        }
                        lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_n(k-Ta,n);
                        //print_DEBUG("F", k+1, " simulated to ", exp(lnF[k]) );
                    } // rates
                } // steps
                // Now have F_k(T_alpha) for k=10,..,M
                F = exp(lnF);
                T floating_swap_rate;
                
                floating_swap_rate = SR_from_F(F, yearly_payments, int_ExTime, int_Tb );
                swap_val = C_ab( F, yearly_payments, int_ExTime, int_Tb) * max(floating_swap_rate - r_fix, 0.0);
                SR[i][t] = floating_swap_rate;
                Libor[i][t] = exp(lnF[ int(exTimes[t]) ]);
                swap_vals[i][t] = swap_val;
            } // paths-loop
        } // Exercise dates-loop
        print("Presim forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
        begin_time = clock();
        
        // Find regression coefficients in backwards loop
        // Exercise value at last exercise time
        vec Y(nPaths_presim);
        vec payoff(nPaths_presim);
        for(int i = 0; i<nPaths_presim; ++i){
            payoff(i) = swap_vals[i][exTimes.size() - 1];
        }
        
        // Backwards loop:
        for(int t = int(exTimes.size() - 2); t >= 1; --t ){
            //print("exTime is ", exTimes[t], " t is ", t);
            vector<double> ITMswapRate, ITMLibor, ITMY;
            vector<int> indices;
            for(int i = 0; i<nPaths_presim; ++i){
                payoff(i) /=  (1 + Libor[i][t]);
            }
            // Find in the money paths
            for(int i = 0; i<nPaths_presim; ++i){
                if(swap_vals[i][t] > 0) { // If ITM
                    ITMswapRate.push_back(SR[i][t]);
                    ITMLibor.push_back(Libor[i][t]);
                    ITMY.push_back(payoff(i));
                    indices.push_back(i);
                };
            }
            size_t ITMsize = ITMswapRate.size();
            mat X(ITMsize, no_basis_funcs, fill::zeros);

            for(int i = 0; i<ITMsize; ++i){
                X.row(i) = vec({ 1.0, ITMswapRate[i], ITMLibor[i] }).t();
            }
            // Do SVD:
            mat U, V;
            vec s;
            svd(U,s,V,X);
            
            mat D(ITMsize, no_basis_funcs, fill::zeros);
            mat Sig(ITMsize, ITMsize, fill::zeros);

            for(int i=0; i<s.n_rows; ++i){
                D(i,i) = s(i);
                Sig(i,i) = 1.0/(s(i)*s(i) + lambda*lambda );
            }
            
            beta.col(t) = V * D.t() * Sig * U.t() * vec(ITMY) ;
            //print("beta col ", t, " is ", beta.col(t));
            
        } // Backward loop
    } // Scope
    print("Presim backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    // Main-simulation
    vector<vector<double>> SR2(nPaths, vector<double>(exTimes.size()));
    vector<vector<double>> Libor2(nPaths, vector<double>(exTimes.size()));
    vector<vector<double>> swap_vals2(nPaths, vector<double>(exTimes.size()));

    // Main-simulation
    for(int t = 1; t<exTimes.size(); ++t ){
        //print("exTime is ", exTimes[t], " t is ", t);
        T swap_val(0.0);
        int int_ExTime = (int)(exTimes[t]);
        int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
        T dt( 1.0/double(nSteps_y) );
        T sqDt = sqrt(dt);
        
        for(int i = 0; i<nPaths; ++i){
            mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
            //print(size(gauss_n));
            vector<T> lnF = log(initF);
            for(int n = 0; n<nSteps; ++n) {
                for(int k = int_ExTime; k < int_Tb; ++k)
                {
                    T sum(0.0);
                    for(int j = int(exTimes[t]); j <= k; ++j ) {
                        double Fj = exp(lnF[j]);
                        sum += corr[k][j] * vol[j][j] * Fj / (1.0 + tau * Fj);
                    }
                    lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_n(k-Ta,n);
                    //print("F", k+1, " simulated to ", exp(lnF[k]) );
                } // rates
            } // steps
            // Now have F_k(T_alpha) for k=10,..,M
            F = exp(lnF);
            T floating_swap_rate;
            floating_swap_rate = SR_from_F(F, yearly_payments, int_ExTime, int_Tb );
            // Compute value of swap
            swap_val = C_ab( F, yearly_payments, int_ExTime, int_Tb) * max(floating_swap_rate - r_fix, 0.0);
            
            SR2[i][t] = floating_swap_rate;
            Libor2[i][t] = exp(lnF[ int_ExTime ]);
            swap_vals2[i][t] = swap_val;
            //print_DEBUG("Set libor to ", Libor2[i][t], " SR ", SR2[i][t], " swap val ",  swap_val, " at time ", t);
        } // paths-loop
    } // Exercise dates-loop
    print("Main forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    
    // Exercise value at last exercise time (exercise if ITM)
    vec payoff_main(nPaths);
    for(int i = 0; i<nPaths; ++i){
        payoff_main(i) = swap_vals2[i][exTimes.size() - 1];
    }
    //print(payoff_main);

    // Backwards loop, from second last exercise time and back :
    for(int t = int(exTimes.size() - 2); t >= 1; --t ){
        vector<double> ITMswapRate, ITMLibor, ITMY;
        vector<int> indices;
        for(int i = 0; i<nPaths; ++i){
            payoff_main(i) /=  (1.0 + Libor2[i][t]);
        }
        // Find in the money paths
        for(int i = 0; i<nPaths; ++i){
            if(swap_vals2[i][t] > 0) { // If ITM
                ITMswapRate.push_back(SR2[i][t]); // Lav evt til vector of vectors
                ITMLibor.push_back(Libor2[i][t]);
                ITMY.push_back(payoff_main(i));
                indices.push_back(i);
            };
        }
        size_t ITMsize = ITMswapRate.size();
        //print("ITMsize is ", ITMsize);
        mat X(ITMsize, no_basis_funcs, fill::zeros);

        for(int i = 0; i<ITMsize; ++i){
            X.row(i) = vec( { 1.0, ITMswapRate[i], ITMLibor[i] }).t();
        }
        vec EY = X*beta.col(t);
        print_DEBUG("X is: ");
        print_DEBUG(X);
        print_DEBUG("EY is ");
        print_DEBUG(EY);
        print_DEBUG("ITMY is ");
        print_DEBUG(ITMY);
        
        for(int i = 0; i<ITMsize; ++i){
            if( swap_vals2[indices[i]][t] > EY[i] ){
                payoff_main(indices[i]) = swap_vals2[indices[i]][t] ;
                print_DEBUG("payoff set to ", payoff_main(indices[i]), " at index ", i, " ", indices[i]);
            }
        }
    } // Exercise times
    print("Main backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
 
    double sumPayoff = sum(payoff_main);
    T disc;
    int int_t = int(t), int_ExT1 = int(exTimes[1]);
    disc = DF_from_F(initF, yearly_payments, int_t, int_ExT1);
    //print("disc is ", disc, " from time ", int(t), " to time ", int(exTimes[1]) );
    return disc * sumPayoff/double(nPaths) * notional;
}



template<class T> T LMM_BermudaSwaptionAAD(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF, vector<double> & exTimes,
                             double t, double Ta, double Tb, T r_fix, double notional,
                             int seed1, int seed2, int nPaths, int nPaths_presim, int nSteps_y, double yearly_payments, int dim_n)
{
    clock_t begin_time = clock();
    T tau(1.0/yearly_payments );
    int int_Ta = (int)(Ta), int_Tb = (int)(Tb);
    
    int no_basis_funcs = 3;
    mat beta(no_basis_funcs, exTimes.size() - 1);
    double lambda = 0.0; // Tikhonov parameter
    
    size_t M1 = initF.size();
    
    vec means(M1, fill::zeros);
    mat C(M1, M1, fill::zeros);
    
    // Fill covariance matrix
    for(int i = 0; i<M1; ++i){
        for(int j = 0; j<M1; ++j){
            C(i,j) = double(corr[i][j]);
        }
    }
    mat corr2 = C( span(int_Ta,M1-1), span(int_Ta,M1-1) );
    // Find eigenvalues and vectors of correlation matrix
    mat P;
    vec eigval;
    eig_sym(eigval, P, corr2);
    // H has diagonal with eigenvalues
    mat H = diagmat(eigval);
    vec largest_n_eigvals = eigval.tail_rows(dim_n);
    // H2 is H but only largest n eigenvalues
    mat H_n = diagmat(largest_n_eigvals);
    // Lambda_n = sqrt(H_n)
    mat L_n = sqrt(H_n);
    // P_n is columns corresponding to largest n eigenvalue
    mat P_n = P.tail_cols(dim_n);
    
    // B_n is used to simulate from
    mat B_n = P_n * L_n;
    
    arma_rng::set_seed(seed1);
    
    vector<double> Fdouble;
    //mrg32k3a myRNG(seed1, seed2) ;
    //myRNG.init(dim_n);
    print("Init took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    
    vector<double> initFdouble;
    for(auto & x : initF ) {initFdouble.push_back(x.value());}
    
    //vector<double> gauss(dim_n);
    { // Scope to destruct everything declared in pre-simulation and first backward loop
        vector<vector<double>> SR(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> Libor(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> swap_vals(nPaths_presim, vector<double>(exTimes.size()));
        
        // Simulate swap rates, pre-simulation
        for(int t = 1; t<exTimes.size(); ++t ){
            //print("exTime is ", exTimes[t], " t is ", t);
            double swap_val(0.0);
            int int_ExTime = (int)(exTimes[t]);
            int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
            double dt( 1.0 /double(nSteps_y));
            double sqDt = sqrt(dt);

            for(int i = 0; i<nPaths_presim; ++i){
                mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
                vector<double> lnF = log(initFdouble);
                for(int n = 0; n<nSteps; ++n) {
                    //print_DEBUG("Time is ", exTimes[t-1] + dt*(n+1));
                    for(int k = int_ExTime; k < int_Tb; ++k)
                    {
                        double sum(0.0);
                        for(int j = int_ExTime; j <= k; ++j ) {
                            double Fj = exp(lnF[j]);
                            sum += corr[k][j].value() * vol[j][j].value() * Fj / (1.0 + tau.value() * Fj);
                        }
                        lnF[k] += vol[k][k].value()*tau.value()*sum*dt - vol[k][k].value()*vol[k][k].value()/2.0*dt + sqDt*vol[k][k].value()*gauss_n(k-Ta,n);
                        //print_DEBUG("F", k+1, " simulated to ", exp(lnF[k]) );
                    } // rates
                } // steps
                // Now have F_k(T_alpha) for k=10,..,M
                Fdouble = exp(lnF);
                double floating_swap_rate;
                
                floating_swap_rate = SR_from_F(Fdouble, yearly_payments, int_ExTime, int_Tb );
                swap_val = C_ab( Fdouble, yearly_payments, int_ExTime, int_Tb) * max(floating_swap_rate - r_fix.value(), 0.0);
                SR[i][t] = floating_swap_rate;
                Libor[i][t] = exp(lnF[ int(exTimes[t]) ]);
                swap_vals[i][t] = swap_val;
            } // paths-loop
        } // Exercise dates-loop
        print("Presim forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
        begin_time = clock();
        
        // Find regression coefficients in backwards loop
        // Exercise value at last exercise time
        vec Y(nPaths_presim);
        vec payoff(nPaths_presim);
        for(int i = 0; i<nPaths_presim; ++i){
            payoff(i) = swap_vals[i][exTimes.size() - 1];
        }
        
        // Backwards loop:
        for(int t = int(exTimes.size() - 2); t >= 1; --t ){
            //print("exTime is ", exTimes[t], " t is ", t);
            vector<double> ITMswapRate, ITMLibor, ITMY;
            vector<int> indices;
            for(int i = 0; i<nPaths_presim; ++i){
                payoff(i) /=  (1 + Libor[i][t]);
            }
            // Find in the money paths
            for(int i = 0; i<nPaths_presim; ++i){
                if(swap_vals[i][t] > 0) { // If ITM
                    ITMswapRate.push_back(SR[i][t]);
                    ITMLibor.push_back(Libor[i][t]);
                    ITMY.push_back(payoff(i));
                    indices.push_back(i);
                };
            }
            size_t ITMsize = ITMswapRate.size();
            mat X(ITMsize, no_basis_funcs, fill::zeros);

            for(int i = 0; i<ITMsize; ++i){
                X.row(i) = vec({ 1.0, ITMswapRate[i], ITMLibor[i] }).t();
            }
            // Do SVD:
            mat U, V;
            vec s;
            svd(U,s,V,X);
            
            mat D(ITMsize, no_basis_funcs, fill::zeros);
            mat Sig(ITMsize, ITMsize, fill::zeros);

            for(int i=0; i<s.n_rows; ++i){
                D(i,i) = s(i);
                Sig(i,i) = 1.0/(s(i)*s(i) + lambda*lambda );
            }
            
            beta.col(t) = V * D.t() * Sig * U.t() * vec(ITMY) ;
            //print("beta col ", t, " is ", beta.col(t));
            
        } // Backward loop
    } // Scope
    print("Presim backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    
    
    
    //number::tape->rewindToMark();
    
    // Main-simulation
    vector<T> F;
    vector<vector<T>> SR2(nPaths, vector<T>(exTimes.size()));
    vector<vector<T>> Libor2(nPaths, vector<T>(exTimes.size()));
    vector<vector<T>> swap_vals2(nPaths, vector<T>(exTimes.size()));
    number::tape->mark();
    // Main-simulation
    for(int t = 1; t<exTimes.size(); ++t ){
        //print("exTime is ", exTimes[t], " t is ", t);
        T swap_val(0.0);
        int int_ExTime = (int)(exTimes[t]);
        int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
        T dt( 1.0/double(nSteps_y) );
        T sqDt = sqrt(dt);
        
        for(int i = 0; i<nPaths; ++i){
            mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
            //print(size(gauss_n));
            vector<T> lnF = log(initF);
            for(int n = 0; n<nSteps; ++n) {
                for(int k = int_ExTime; k < int_Tb; ++k)
                {
                    T sum(0.0);
                    for(int j = int(exTimes[t]); j <= k; ++j ) {
                        T Fj(exp(lnF[j]));
                        sum += corr[k][j] * vol[j][j] * Fj / (1.0 + tau * Fj);
                    }
                    lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_n(k-Ta,n);
                    //print("F", k+1, " simulated to ", exp(lnF[k]) );
                } // rates
            } // steps
            // Now have F_k(T_alpha) for k=10,..,M
            F = exp(lnF);
            T floating_swap_rate;
            floating_swap_rate = SR_from_F(F, yearly_payments, int_ExTime, int_Tb );
            // Compute value of swap
            swap_val = C_ab( F, yearly_payments, int_ExTime, int_Tb) * max(floating_swap_rate - r_fix, 0.0);
            
            SR2[i][t] = floating_swap_rate;
            Libor2[i][t] = exp(lnF[ int_ExTime ]);
            swap_vals2[i][t] = swap_val;
            //swap_vals2[i][t].propagateToMark();
            //print_DEBUG("Set libor to ", Libor2[i][t], " SR ", SR2[i][t], " swap val ",  swap_val, " at time ", t);
        } // paths-loop
    } // Exercise dates-loop
    print("Main forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    
    // Exercise value at last exercise time (exercise if ITM)
    //vec payoff_main(nPaths);
    vector<T> payoff_main(nPaths);
    for(int i = 0; i<nPaths; ++i){
        payoff_main[i] = swap_vals2[i][exTimes.size() - 1];
    }

    // Backwards loop, from second last exercise time and back :
    for(int t = int(exTimes.size() - 2); t >= 1; --t ){
        vector<T> ITMswapRate, ITMLibor, ITMY;
        vector<int> indices;
        for(int i = 0; i<nPaths; ++i){
            payoff_main[i] /=  (1.0 + Libor2[i][t]);
        }
        // Find in the money paths
        for(int i = 0; i<nPaths; ++i){
            if(swap_vals2[i][t] > 0) { // If ITM
                ITMswapRate.push_back(SR2[i][t]); // Lav evt til vector of vectors
                ITMLibor.push_back(Libor2[i][t]);
                ITMY.push_back(payoff_main[i]);
                indices.push_back(i);
            };
        }
        size_t ITMsize = ITMswapRate.size();
        //print("ITMsize is ", ITMsize);
        mat X(ITMsize, no_basis_funcs, fill::zeros);

        for(int i = 0; i<ITMsize; ++i){
            X.row(i) = vec( { 1.0, ITMswapRate[i].value(), ITMLibor[i].value() }).t();
        }
        vec EY = X*beta.col(t);
        /*print_DEBUG("X is: ");
        print_DEBUG(X);
        print_DEBUG("EY is ");
        print_DEBUG(EY);
        print_DEBUG("ITMY is ");
        print_DEBUG(ITMY);*/
        
        for(int i = 0; i<ITMsize; ++i){
            if( swap_vals2[indices[i]][t] > EY[i] ){
                payoff_main[indices[i]] = swap_vals2[indices[i]][t] ;
                //print_DEBUG("payoff set to ", payoff_main[indices[i]], " at index ", i, " ", indices[i]);
            }
        }
    } // Exercise times
    print("Main backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
 
    T sumPayoff( sum(payoff_main) );
    T disc;
    int int_t = int(t), int_ExT1 = int(exTimes[1]);
    disc = DF_from_F(initF, yearly_payments, int_t, int_ExT1);
    //print("disc is ", disc, " from time ", int(t), " to time ", int(exTimes[1]) );
    T res( disc * sumPayoff/double(nPaths) * notional ); 
    //number::propagateMarkToStart();
    return res;
}


#endif /* Bermuda_h */
