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
#include "/usr/local/Cellar/armadillo/9.700.2/include/armadillo"

using namespace std;

template<class T> T LMM_swaption(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF,
                             T t, T Ta, T Tb, T r_fix, T notional,
                             int seed1, int seed2, int nPaths, int nSteps, T yearly_payments)
{
    
    T dt(double(Ta-t)/double(nSteps));
    T tau(1.0/yearly_payments );
    T sqDt = sqrt(dt);
    T res(0.0);
    T avg_float(0.0);
    size_t M = initF.size();
    
    T nPayments( (Tb - Ta) * yearly_payments);
    print("payments ", nPayments);
    
    vec means(M, fill::zeros);
    mat C(M, M, fill::zeros);
    
    // Fill covariance matrix
    for(int i = 0; i<M; ++i){
        for(int j = 0; j<M; ++j){
            //C(i,j) = corr[i][j] * sqrt(vol[i][i] * vol[j][j]);
            C(i,j) = corr[i][j] * vol[i][i] * vol[j][j];
            //print(i, " ", j, " Corr ", corr[i][j], " vol i ", sqrt(vol[i][i]), " vol j ", sqrt(vol[j][j]), " gives ", C(i,j) );
        }
    }
    //print("Cov is");
    //print(C);
    
    arma_rng::set_seed(seed1);
    
    vector<T> F;
    int count = 0;
    
    for(int i = 0; i<nPaths; ++i){
        // Simulate multivariate gaussians with covariance C and mean 0
        mat gauss = arma::mvnrnd(means, C, nSteps);
        
        vector<T> lnF = log(initF);
        
        for(int n = 0; n<nSteps; ++n) {
            //print("STEP ", n);
            for(int k = int(Ta); k<M; ++k) // Loop over number of Forward rates to simulate, 0,1,2 = F1,F2,F3
            {
                //print("k is ", k);
                // Compute sum in Brigo 6.53
                T sum(0.0);
                for(int j = 0; j <= k; ++j ) { // Loop over yearly forward rates,
                    // 2.0 = alpha+1 = k0, 3.0 = alpha+2 = k1, 4.0 = alpha+3 = beta = k2,
                    //print("j is ", j);
                    double Fj = exp(lnF[j]);
                    sum += corr[k][j] * tau * vol[j][j] * Fj / (1.0 + tau * Fj);
                }
                //print("sum is ", sum);
                lnF[k] += vol[k][k] * sum * dt - vol[k][k] * vol[k][k]/2.0 * dt + sqDt * gauss(k,n); // vol[k][k] *
                //print("F[k] is ", exp(lnF[k]));
            }
            //print(exp(lnF));
        }
        // Now have F_k(T_alpha) for k=1,..,M
        F = exp(lnF);
        //print("Sim F is ");
        //print(F);
        T floating_swap_rate;
        floating_swap_rate = SR_from_F(F, yearly_payments, (int)(Ta), int(Tb) );
        //print("floating_swap_rate is ", floating_swap_rate);
        
        T diff_rates;
        diff_rates = max(floating_swap_rate - r_fix, 0.0);  // Payer: floating_swap_rate - r_fix
        count += diff_rates > 0.0 ? 1 : 0;
        
        // Compute value of swap using this floating rate
        T val(0.0);
        for(int j=0; j<nPayments; ++j){
            // Discount back to time Ta
            T disc = DF_from_F(F, yearly_payments, int(Ta), int(Ta) + 1 + j);
            //print("Disc is ", disc);
            val += disc * diff_rates * notional;
        }
        //print("val is ", val);
        res += val;

        avg_float += floating_swap_rate;
    }
    print("avg_float is ", avg_float/double(nPaths));
    print("count is ", count);
    // Discounting from Ta to t:
    T disc;
    disc = DF_from_F(initF, yearly_payments, int(t), int(Ta));
    print("Disc is ", disc);
    return disc * res/double(nPaths);
}


#endif /* LMM_h */

/*
 T lnF1 = log(initF[0]);
 T lnF2 = log(initF[1]);
 T lnF3 = log(initF[2]);
 lnF1 += vol[0][0] *
corr[0][0] * 1/yearly_payments * vol[0][0] * exp(lnF1) / (1 + 1/yearly_payments * exp(lnF1)) * dt - vol[0][0] * vol[0][0]/2.0 * dt + vol[0][0] * sqDt * gauss[n];
                   
lnF2 += vol[1][1] *
( corr[1][0] * 1/yearly_payments * vol[0][0] * exp(lnF1) / (1 + 1/yearly_payments * exp(lnF1))
 +
  corr[1][1] * 1/yearly_payments * vol[1][1] * exp(lnF2) / (1 + 1/yearly_payments * exp(lnF2))
 ) * dt
- vol[1][1] * vol[1][1]/2.0 * dt + vol[1][1] * sqDt * gauss[n + nSteps];

lnF3 += vol[2][2] *
( corr[2][0] * 1/yearly_payments * vol[0][0] * exp(lnF1) / (1 + 1/yearly_payments * exp(lnF1))
 +
  corr[2][1] * 1/yearly_payments * vol[1][1] * exp(lnF2) / (1 + 1/yearly_payments * exp(lnF2))
 +
  corr[2][2] * 1/yearly_payments * vol[2][2] * exp(lnF3) / (1 + 1/yearly_payments * exp(lnF3))
 ) * dt
- vol[2][2] * vol[2][2]/2.0 * dt + vol[2][2] * sqDt * gauss[n + nSteps*2];
cout << exp(lnF1) << " " << exp(lnF2) << " " << exp(lnF3) << endl;*/
