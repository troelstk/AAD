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
    size_t M = int_Tb-int_Ta;
    
    vector<vector<T>> cov_s(M, vector<T>(M));
    vector<vector<T>> lower(M, vector<T>(dim_n));
    vector<vector<double>> cov_d(M, vector<double>(M));
    vector<vector<double>> lower_d(M, vector<double>(dim_n));
    
    mat Cov(M, M, fill::zeros);
    // Fill covariance vec of vecs
    for(int i = int_Ta; i<Tb; ++i){
        for(int j = int_Ta; j<Tb; ++j){
            T cov_ij(corr[i][j]);
            cov_s[i-int_Ta][j-int_Ta] = cov_ij;
            cov_d[i-int_Ta][j-int_Ta] = double(cov_ij);
            Cov(i-int_Ta, j-int_Ta) = double(cov_ij);
        }
    }
    mat P;
    vec eigval;
    eig_sym(eigval, P, Cov);
    vec largest_n_eigvals = eigval.tail_rows(dim_n);
    mat H_n = diagmat(largest_n_eigvals);
    mat L_n = sqrt(H_n);
    mat P_n = P.tail_cols(dim_n);
    mat B_n = P_n * L_n;
    
 
    
    vector<vector<T>> Pert =
    {{T(1.0), T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0)},
    {T(0.0), T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(1.0)},
    {T(0.0), T(0.0),T(1.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0)},
    {T(0.0), T(0.0),T(0.0),T(0.0),T(0.0),T(1.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0)},
    {T(0.0), T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(1.0),T(0.0),T(0.0),T(0.0),T(0.0)},
    {T(0.0), T(0.0),T(0.0),T(0.0),T(1.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0)},
    {T(0.0), T(1.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0)},
    {T(0.0), T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(1.0),T(0.0),T(0.0),T(0.0)},
    {T(0.0), T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(1.0),T(0.0),T(0.0)},
    {T(0.0), T(0.0),T(0.0),T(1.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0)},
    {T(0.0), T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(0.0),T(1.0),T(0.0)}};
    
    vector<vector<T>> L_num(dim_n, vector<T>(dim_n, T(0.0)));
    vector<vector<T>> eig_num(dim_n, vector<T>(dim_n, T(0.0)));
    vector<vector<T>> P_num(M, vector<T>(dim_n, T(0.0))); // (rows, vector<T>(cols, T(0.0)));
    vector<vector<T>> B_num(M, vector<T>(dim_n, T(0.0)));

    
    for(int i = 0; i<dim_n; ++i){
        eig_num[i][i] = T(largest_n_eigvals(i));
        L_num[i][i] = sqrt(eig_num[i][i]);
    }
    
    // Copy eigen vectors of largest eigenvalues to overloaded P
    for(int i = 0; i<M; i++){ // Rows
        for(int j = 0; j<dim_n; j++){ // cols
            P_num[i][j] = T(P_n(i,j));
        }
    }

    B_num = MatMatProd(P_num, L_num);
    
    vector<vector<T>> A = MatMatProd(MatMatProd(transp(Pert), cov_s), Pert);
    
    lower = CholMod(A, B_num);
    

    for(int i = 0; i<M; i++){ // Rows
        for(int j = 0; j<dim_n; j++){ // cols
            lower_d[i][j] = double(lower[i][j]);
        }
    }
    
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
    for(auto & x : initF ) {initFdouble.push_back(double(x));}
    
    { // Scope to destruct everything declared in pre-simulation and backward loop
        mrg32k3a myRNG(seed1, seed2);
        myRNG.init(dim_n);
        vector<double> gauss(dim_n);
        vector<vector<double>> SR(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> Libor(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> swap_vals(nPaths_presim, vector<double>(exTimes.size()));
        
        // Simulate swap rates, pre-simulation
        double dt( 1.0 /double(nSteps_y));
        double sqDt = sqrt(dt);

        for(int i = 0; i<nPaths_presim; ++i) {
            for(int t = 1; t<exTimes.size(); ++t ) {
                double swap_val(0.0);
                int int_ExTime = (int)(exTimes[t]);
                int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
                vector<double> lnF = log(initFdouble);
                
                for(int n = 0; n<nSteps; ++n) {
                    myRNG.nextG(gauss);
                    vector<double> gauss_corr = MatVecProd( lower_d, gauss);
                    for(int k = int_ExTime; k < int_Tb; ++k)
                    {
                        double sum(0.0);
                        for(int j = int_ExTime; j <= k; ++j ) {
                            double Fj = exp(lnF[j]);
                            sum += double(corr_vol[k][j]) * Fj / (1.0 + double(tau) * Fj);
                        }
                        lnF[k] += double(vol[k][k])*double(tau)*sum*dt - double(vol[k][k])*double(vol[k][k])/2.0*dt + sqDt*double(vol[k][k])*gauss_corr[k-int_Ta];
                    } // rates
                } // steps
                // Now have F_k(T_alpha) for k=10,..,M
                Fdouble = exp(lnF);
                double floating_swap_rate = SR_from_F(Fdouble, yearly_payments, int_ExTime, int_Tb );
                
                swap_val = double(disc) * notional *
                    C_ab( Fdouble, yearly_payments, int_ExTime, int_ExTime, int_Tb) *
                    max(floating_swap_rate - double(r_fix), 0.0);
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
        for(int t = int(exTimes.size() - 2); t >= 1; --t )
        {
            vector<double> ITMswapRate, ITMLibor, ITMY;
            vector<int> indices;
            for(int i = 0; i<nPaths_presim; ++i){
                payoff(i) /=  (1.0 + Libor[i][t]);
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
            mat EY = X*beta.col(t);
            for(int i = 0; i<indices.size(); ++i){
                if(swap_vals[indices[i]][t] > EY(i) ) { // If ITM
                    payoff(i) = swap_vals[indices[i]][t];
                };
            }
        } // Backward loop
    } // Scope
    print_DEBUG("Presim backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    
    // Main-simulation
    vector<T> F;

    T dt( 1.0/double(nSteps_y) );
    T sqDt = sqrt(dt);

    // Main-simulation
    T end_result(0.0);
    
    vector<double> swap_rates(nPaths);
    mrg32k3a myRNG(seed2, seed1);
    myRNG.init(dim_n);
    vector<double> gauss(dim_n);

    vector<int> count_ex = vector<int>(exTimes.size(), 0);
    vector<double> ex_boundary = vector<double>(exTimes.size(), 100);
    
    
    number::tape->mark();
    for(int i = 0; i<nPaths; ++i){
        number::tape->rewindToMark();
        int eta = 1;
        T swap_path(0.0);
        for(int t = 1; t<exTimes.size(); ++t ){
            int int_ExTime = (int)(exTimes[t]);
            int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
            T swap_val(0.0);
            vector<T> lnF = log(initF);
            for(int n = 0; n<nSteps; ++n) {
                myRNG.nextG(gauss);
                vector<T> gauss_corr = MatVecProd(lower, gauss);
                for(int k = int_ExTime; k < int_Tb; ++k)
                {
                    T sum(0.0);
                    for(int j = int_ExTime; j <= k; ++j ) {
                        T Fj(exp(lnF[j]));
                        sum += corr_vol[k][j] * Fj / (1.0 + tau * Fj);
                    }
                    lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_corr[k-int_Ta];
                } // rates
            } // steps
            // Now have F_k(T_alpha) for k=10,..,M
            F = exp(lnF);

            T floating_swap_rate;
            floating_swap_rate = SR_from_F(F, yearly_payments, int_ExTime, int_Tb );

            swap_val = disc * notional *
                C_ab( F, yearly_payments, int_ExTime, int_ExTime, int_Tb) *
                max(floating_swap_rate - r_fix, 0.0);

            swap_path += double(eta) * swap_val;
            
            double EY2 = t < exTimes.size() - 1 ?
                // continuation value - at last ex time the eta value computed below will never be used (so 0.0 could be anything), just avoid beta of last t
                as_scalar(vec( { 1.0, double(floating_swap_rate), double(F[ int_ExTime ]) }).t() * beta.col(t)) : 0.0;

            eta *= EY2 > swap_val ? 1 : 0; // 1 if continue, 0 if exercise
        }
        swap_path.propagateToMark();
        
        end_result += swap_path;
    }
    print_DEBUG("Main forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");

   
    number::propagateMarkToStart();
    
    
    return end_result/double(nPaths) ;
}

#endif /* BermudaAAD_h */
