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
#include "Cholesky.h"

template<class T> T LMM_BermudaSwaptionAAD(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF, vector<double> & exTimes,
                             double t, double Ta, double Tb, T r_fix, double notional,
                             int seed1, int seed2, int nPaths, int nPaths_presim, int nSteps_y, double yearly_payments, int dim_n)
{
    T tau(1.0/yearly_payments );
    int int_Ta = (int)(Ta), int_Tb = (int)(Tb);
    int no_basis_funcs = 3;
    mat beta(no_basis_funcs, exTimes.size() - 1);
    double lambda = 0.0; // Tikhonov parameter
    
    size_t M1 = initF.size();
    
    double eps = 0.01;

    mat C(M1, M1, fill::zeros);
    mat Cov(M1, M1, fill::zeros);
    print("Cov is ");
    // Fill covariance matrix
    for(int i = 0; i<M1; ++i){
        for(int j = 0; j<M1; ++j){
            C(i,j) = double(corr[i][j]);
            double entry = double(corr[i][j]) * double(vol[i][i]) * double(vol[j][j]);
            Cov(i,j) = i == j ? entry + eps : entry ;
            //cout << C(i,j) << ",";
        }
        //cout << "\n";
    }
    
    vector<vector<T>> cov_s(int_Tb-int_Ta, vector<T>(int_Tb-int_Ta));
    vector<vector<T>> lower(int_Tb-int_Ta, vector<T>(int_Tb-int_Ta));
    
    // test: Fill covariance vec of vecs
    for(int i = int_Ta; i<Tb; ++i){
        for(int j = int_Ta; j<Tb; ++j){
            T cov_ij(corr[i][j] * double(vol[i][i]) * double(vol[j][j]));
            cov_s[i-int_Ta][j-int_Ta] = (i == j) ? cov_ij + eps : cov_ij ;
            
            //cout << cov_s[i-int_Ta][j-int_Ta].value() << ",";
            //cout << i << " " << j << "\n";
        }
        //cout <<"\n" ;
    }
    
    lower = myChol(cov_s);
    
    
    mat Cov_small = Cov( span(int_Ta,M1-1), span(int_Ta,M1-1) );
    print(Cov_small);
    
    print(Cov_small.is_symmetric(1e-5));
    print(Cov_small.is_sympd(1e-5));
    
    print(arma::chol(Cov_small));
    
    //print(size(corr2));
    mat corr2 = C( span(int_Ta,M1-1), span(int_Ta,M1-1) );
    // Find eigenvalues and vectors of correlation matrix
    mat P;
    vec eigval;
    eig_sym(eigval, P, corr2);
    // H has diagonal with eigenvalues
    //mat H = diagmat(eigval);
    vec largest_n_eigvals = eigval.tail_rows(dim_n);
    // H_n is diagonal with only largest n eigenvalues
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
    mrg32k3a myRNG(seed1, seed2);
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
        double dt( 1.0 /double(nSteps_y));
        double sqDt = sqrt(dt);

        for(int i = 0; i<nPaths_presim; ++i) {
            for(int t = 1; t<exTimes.size(); ++t ) {
                //print("exTime is ", exTimes[t], " t is ", t);
                double swap_val(0.0);
                int int_ExTime = (int)(exTimes[t]);
                int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
                
                /*myRNG.init(nSteps);
                vector<double> gauss(nSteps);
                
                myRNG.nextG(gauss);*/
                
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
                        //lnF[k] += vol[k][k].value()*tau.value()*sum*dt - vol[k][k].value()*vol[k][k].value()/2.0*dt + sqDt*vol[k][k].value()*gauss[n];;
                        //print_DEBUG("F", k+1, " simulated to ", exp(lnF[k]) );
                    } // rates
                } // steps
                // Now have F_k(T_alpha) for k=10,..,M
                Fdouble = exp(lnF);
                double floating_swap_rate = SR_from_F(Fdouble, yearly_payments, int_ExTime, int_Tb );
                
                swap_val = disc.value() * notional *
                    C_ab( Fdouble, yearly_payments, int_ExTime, int_ExTime, int_Tb) *
                    max(floating_swap_rate - r_fix.value(), 0.0);
                SR[i][t] = floating_swap_rate;
                Libor[i][t] = Fdouble[ int_ExTime ];
                swap_vals[i][t] = swap_val;
                
            } // exTimes
        } // paths
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
    
    arma_rng::set_seed(seed2);
    // Main-simulation
    vector<T> F;

    T dt( 1.0/double(nSteps_y) );
    T sqDt = sqrt(dt);

    // Main-simulation
    T end_result(0.0);
    

    number::tape->mark();
    for(int i = 0; i<nPaths; ++i){
        number::tape->rewindToMark();
        double eta = 1;
        T swap_path(0.0);
        for(int t = 1; t<exTimes.size(); ++t ){
            //print("exTime is ", exTimes[t], " t is ", t);
            int int_ExTime = (int)(exTimes[t]);
            int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
            
            /*myRNG.init(nSteps);
            vector<double> gauss(nSteps);
            
            myRNG.nextG(gauss);*/
            T swap_val(0.0);
            mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
            
            //mat gauss_M = arma::mvnrnd(zeros(int_Tb-int_Ta), corr2, nSteps);
            vector<T> lnF = log(initF);
            for(int n = 0; n<nSteps; ++n) {
                for(int k = int_ExTime; k < int_Tb; ++k)
                {
                    T sum(0.0);
                    for(int j = int_ExTime; j <= k; ++j ) {
                        T Fj(exp(lnF[j]));
                        sum += corr_vol[k][j] * Fj / (1.0 + tau * Fj);
                    }
                    //lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss[n];
                    lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_n(k-Ta,n);
                    //lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_M(k-Ta,n);
                } // rates
            } // steps
            // Now have F_k(T_alpha) for k=10,..,M
            F = exp(lnF);
            T floating_swap_rate;
            floating_swap_rate = SR_from_F(F, yearly_payments, int_ExTime, int_Tb );
            
            swap_val = disc * notional *
                C_ab( F, yearly_payments, int_ExTime, int_ExTime, int_Tb) *
                max(floating_swap_rate - r_fix, 0.0);

            swap_path += eta * swap_val;
            
            double EY2 = t < exTimes.size() - 1 ?
                // continuation value - at last ex time the eta value computed below will never be used (so 0.0 could be anything)
                as_scalar(vec( { 1.0, floating_swap_rate.value(), F[ int_ExTime ].value() }).t() * beta.col(t)) : 0.0;
            eta *= EY2 > swap_val ? 1.0 : 0.0; // 1 if continue, 0 if exercise
        }
        swap_path.propagateToMark();
        end_result += swap_path;
    }
    print_DEBUG("Main forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    number::propagateMarkToStart();

    return end_result/double(nPaths) ;
}

#endif /* BermudaAAD_h */
