//
//  LMM.h
//  AAD
//
//  Created by Troels Tang Karlsen on 22/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef LMM_h
#define LMM_h
#include <cmath>
#include "utilities.h"
#include "gaussians.h"
#include "mrg32k3a.h"
#include "/usr/local/Cellar/armadillo/9.700.2/include/armadillo"

using namespace std;

template<class T> T LMM_swaption(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF,
                             T t, T Ta, T Tb, T r_fix, T notional,
                             int seed1, int nPaths, int nSteps, T yearly_payments, int dim_n)
{
    
    T dt(double(Ta-t)/double(nSteps));
    T tau(1.0/yearly_payments );
    T sqDt = sqrt(dt);
    T res(0.0);
    T avg_val(0.0);
    size_t M1 = initF.size();
    
    vec means(M1, fill::zeros);
    mat C(M1, M1, fill::zeros);
    mat C2(M1, M1, fill::zeros);
    
    // Fill covariance matrix
    for(int i = 0; i<M1; ++i){
        for(int j = 0; j<M1; ++j){
            C(i,j) = corr[i][j] * vol[i][i] * vol[j][j];
            C2(i,j) = corr[i][j];
            //print(i, " ", j, " Corr ", corr[i][j], " vol i ", vol[i][i], " vol j ", vol[j][j], " gives ", C(i,j) );
        }
    }

    mat Cs = C( span(int(Ta),M1-1), span(int(Ta),M1-1) ); // Select subset of full covariance matrix to simulate from
    mat corr2 = C2( span(int(Ta),M1-1), span(int(Ta),M1-1) );
    // Find eigenvalues and vectors of correlation matrix
    mat P;
    vec eigval;
    eig_sym(eigval, P, corr2);
    // H has diagonal with eigenvalues
    mat H = diagmat(eigval);
    
    // Choose rank n:
    if( eigval.is_sorted() == 0 ) return 0.0;
    // Assumes eigval is sorted ascending
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
    mrg32k3a myRNG(seed1, 123) ;
    myRNG.init(dim_n);

    vector<double> gauss(dim_n);
    vector<T> lnRates(nPaths);
    
    for(int i = 0; i<nPaths; ++i){
        mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
        vector<T> lnF = log(initF);
        for(int n = 0; n<nSteps; ++n) {
            for(int k = int(Ta); k < int(Tb); ++k) // Loop over number of Forward rates to simulate, 0=F1, 9=F10,..,19=F20
            {   // When Ta is 9, we access 10'th (alpha+1) entry in F
                // Compute sum in Brigo 6.53
                T sum(0.0);
                for(int j = int(Ta); j <= k; ++j ) { // Loop over yearly forward rates,
                    // 2.0 = F_alpha+1 = F_2 = k[1], 3.0 = alpha+2 = k[2], 4.0 = alpha+3 = beta = k[3],
                    double Fj = exp(lnF[j]);
                    sum += corr[k][j] * vol[j][j] * Fj / (1.0 + tau * Fj);
                }
                lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_n(k-Ta,n);
            }
        }
        // Now have F_k(T_alpha) for k=10,..,M
        F = exp(lnF);
        T floating_swap_rate;
        floating_swap_rate = SR_from_F(F, yearly_payments, (int)(Ta), (int)(Tb) );

        lnRates[i] = log(floating_swap_rate);
        
        T diff_rates;
        diff_rates = max(floating_swap_rate - r_fix, 0.0);  // Payer: floating_swap_rate - r_fix
        
        // Compute value of swap using this floating rate
        double C = C_ab( F, yearly_payments, int(Ta), int(Tb));
        res += C * diff_rates;

        avg_val += C;
    }
    double avg = avg_val/double(nPaths);
    print("avg_val is ", avg);

    print("std is ", stdev(lnRates));
    print("kurtosis is ", kurtosis(lnRates));
    print("skew is ", skew(lnRates));
    
    // Discounting back to time t:
    T disc;
    disc = DF_from_F(initF, yearly_payments, int(t), int(Ta));
    return disc * res/double(nPaths) * notional;
}


#endif /* LMM_h */

