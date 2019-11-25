//
//  BermudaAAD.h
//  AAD
//
//  Created by Troels Tang Karlsen on 18/11/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef BermudaAAD_h
#define BermudaAAD_h

#include <cmath>
#include "utilities.h"
#include "gaussians.h"
#include "mrg32k3a.h"
#include "/usr/local/Cellar/armadillo/9.700.2/include/armadillo"


template<class T> T LMM_BermudaSwaptionAAD(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF, vector<double> & exTimes,
                             double t, double Ta, double Tb, T r_fix, double notional,
                             int seed1, int seed2, int nPaths, int nPaths_presim, int nSteps_y, double yearly_payments, int dim_n)
{
    T tau(1.0/yearly_payments );
    int int_Ta = (int)(Ta), int_Tb = (int)(Tb);
    int idx  = 18;
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
    
    // Pre compute product, instead of doing it in inner most loop
    vector<vector<T>> corr_vol(M1, vector<T>(M1));
    for(int k = 1; k < int_Tb; ++k)
    {
        for(int j = 1; j <= k; ++j ) {
            corr_vol[k][j] = corr[k][j] * vol[j][j];
        }
    }
    
    vector<double> Fdouble;

    clock_t begin_time = clock();
    
    T disc;
    int int_t = int(t), int_ExT1 = int(exTimes[1]);
    disc = DF_from_F(initF, yearly_payments, int_t, int_ExT1);
    
    vector<double> initFdouble;
    for(auto & x : initF ) {initFdouble.push_back(x.value());}
    
    { // Scope to destruct everything declared in pre-simulation and first backward loop
        arma_rng::set_seed(seed1);
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
                            sum += corr_vol[k][j].value() * Fj / (1.0 + tau.value() * Fj);
                        }
                        lnF[k] += vol[k][k].value()*tau.value()*sum*dt - vol[k][k].value()*vol[k][k].value()/2.0*dt + sqDt*vol[k][k].value()*gauss_n(k-Ta,n);
                        //print_DEBUG("F", k+1, " simulated to ", exp(lnF[k]) );
                    } // rates
                } // steps
                // Now have F_k(T_alpha) for k=10,..,M
                Fdouble = exp(lnF);
                double floating_swap_rate;
                
                floating_swap_rate = SR_from_F(Fdouble, yearly_payments, int_ExTime, int_Tb );
                swap_val = disc.value() * notional *
                    C_ab( Fdouble, yearly_payments, int_ExTime, int_Tb) *
                    max(floating_swap_rate - r_fix.value(), 0.0);
                SR[i][t] = floating_swap_rate;
                Libor[i][t] = Fdouble[ int_ExTime ];
                swap_vals[i][t] = swap_val;
            } // paths-loop
        } // Exercise dates-loop
        print_DEBUG("Presim forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
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
    print_DEBUG("Presim backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    //print(beta);

    
    arma_rng::set_seed(seed2);
    // Main-simulation
    vector<T> F;
    vector<vector<T>> SR2(nPaths, vector<T>(exTimes.size()));
    vector<vector<T>> Libor2(nPaths, vector<T>(exTimes.size()));
    vector<vector<T>> swap_vals2(nPaths, vector<T>(exTimes.size()));

    T dt( 1.0/double(nSteps_y) );
    T sqDt = sqrt(dt);

    number::tape->mark();
    // Main-simulation
    for(int t = 1; t<exTimes.size(); ++t ){
        //print("exTime is ", exTimes[t], " t is ", t);
        int int_ExTime = (int)(exTimes[t]);
        int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
        for(int i = 0; i<nPaths; ++i){
            number::tape->rewindToMark();
            T swap_val(0.0);
            mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
            //print(size(gauss_n));
            vector<T> lnF = log(initF);
            for(int n = 0; n<nSteps; ++n) {
                for(int k = int_ExTime; k < int_Tb; ++k)
                {
                    T sum(0.0);
                    for(int j = int_ExTime; j <= k; ++j ) {
                        T Fj(exp(lnF[j]));
                        sum += corr_vol[k][j] * Fj / (1.0 + tau * Fj);
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
            swap_val = disc * notional *
                C_ab( F, yearly_payments, int_ExTime, int_Tb) *
                max(floating_swap_rate - r_fix, 0.0);
            SR2[i][t] = floating_swap_rate;
            Libor2[i][t] = F[ int_ExTime ];
            
            //double EY2 = t < exTimes.size() - 1 ?
            //    as_scalar(vec( { 1.0, SR2[i][t].value(), Libor2[i][t].value() }).t() * beta.col(t)) : 0;
            
            swap_val.propagateToMark();
            swap_vals2[i][t] = swap_val;
            //print_DEBUG("Set libor to ", Libor2[i][t], " SR ", SR2[i][t], " swap val ",  swap_val, " at time ", t);
        } // paths-loop
        //print_DEBUG("Adj is ", vol[idx][idx].adjoint());
    } // Exercise dates-loop
    print_DEBUG("Main forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    //print_DEBUG("Adj is ", vol[idx][idx].adjoint());
    //number::propagateMarkToStart();
    
    //vector<double> continues(nPaths);
    //std::fill(continues.begin(), continues.end(), 1.0);
    
    
    // Exercise value at last exercise time (exercise if ITM)
    vector<T> payoff_main(nPaths);
    for(int i = 0; i<nPaths; ++i){
        payoff_main[i] = swap_vals2[i][exTimes.size() - 1];
    }
    double exercises = 0.0;
    //number::tape->mark();
    // Backwards loop, from second last exercise time and back :
    for(int t = int(exTimes.size() - 2); t >= 1; --t ){
        
        for(int i = 0; i<nPaths; ++i){
            //number::tape->rewindToMark();
            payoff_main[i] /=  (1.0 + Libor2[i][t]);
            double EY2 = as_scalar(vec( { 1.0, SR2[i][t].value(), Libor2[i][t].value() }).t() * beta.col(t));
            double exercise = swap_vals2[i][t] > EY2 ? 1.0 : 0.0;
            exercises += exercise;
            T payoff;
            payoff = exercise * swap_vals2[i][t] + (1 - exercise) * payoff_main[i];
            //payoff_main[i] = exercise * swap_vals2[i][t] + (1 - exercise) * payoff_main[i];
            //payoff.propagateToMark();
            payoff_main[i] = payoff;
        }
        //print("Adj is ", vol[idx][idx].adjoint());
        //print("There were ", exercises, " exercises");
    } // Exercise times
    
    
    
    T sumPayoff( sum(payoff_main) );
    //sumPayoff.propagateToMark();
    
    //print("disc is ", disc, " from time ", int(t), " to time ", int(exTimes[1]) );
    T res( sumPayoff/double(nPaths) );

    //number::propagateMarkToStart();
    //res.propagateToStart();
    
    print("disc is ", disc.value());
    print_DEBUG("Main backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    
    return res;
}

#endif /* BermudaAAD_h */
