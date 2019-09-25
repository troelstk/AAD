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
                             int seed1, int seed2, int nPaths, int nSteps, T yearly_payments, T strike)
{
    
    T dt((Ta-t)/double(nSteps));
    T sqDt = sqrt(dt);
    T res(0.0);
    T avg_float(0.0);
    int M = initF.size();
    
    T nPayments( (Tb - Ta) * yearly_payments + 1);
    //print("payments ", nPayments);
    
    // Initialize RNG to M (initF.size()), so that we can simulate all forward rates using same RNG
    mrg32k3a myRNG(seed1, seed2);
    T nGauss(nSteps*initF.size());
    myRNG.init(nGauss);
    vector<double> gauss(nGauss);
    
    vec means(initF.size(), fill::zeros);
    mat C(M, M, fill::zeros);
    
    // Fill cov
    for(int i = 0; i<M; ++i){
        for(int j = 0; j<M; ++j){
            C(i,j) = corr[i][j]*sqrt(vol[i][i])*sqrt(vol[j][j]);
        }
    }
    //print("Cov is ", C);
    arma_rng::set_seed(seed1);
    mat X = arma::mvnrnd(means, C, 5);
    
    print("Simulated values is \n", X);
    
    vector<T> F;
    
    for(int i = 0; i<nPaths; ++i){
        myRNG.nextG(gauss);
        vector<T> lnF = log(initF);
        
        for(int n = 0; n<nSteps; ++n) {
            for(int k = 0; k<Tb-Ta; ++k) // Loop over number of Forward rates to simulate, 0,1,2 = F1,F2,F3
            {
                T sum(0.0);
                for(int j = 0; j <= k; ++j ) { // Loop over yearly forward rates
                    sum += corr[k][j] * 1.0/yearly_payments * vol[j][j] * exp(lnF[j]) /
                    (1.0 + 1.0/yearly_payments * exp(lnF[j]));
                }
                lnF[k] += vol[k][k] * sum * dt - vol[k][k] * vol[k][k]/2.0 * dt +
                vol[k][k] * sqDt *  gauss[n + k * nSteps]; // invNormalCdf(double(rand() + 1.0) / double(RAND_MAX + 2.0));
            }

        }
        // Now have F_k(T_alpha) for k=
        F = exp(lnF);
        T floating_swap_rate = SR_from_F(F, yearly_payments, Ta, Tb);
        
        T val = 0.0;
        
        for(int j=0; j<nPayments; ++j){
            T Ti(Ta + double(j)/yearly_payments);
            //print("Ti is ", Ti);
            
            T disc = DF_from_F(F, yearly_payments, Ta, Ti);
            //print("disc is ", disc);
            
            val += disc * (r_fix - floating_swap_rate) * notional;
        }
        //print("val is ", val);
        res += max(val-strike, 0.0);
        print("float is ", floating_swap_rate);
        avg_float += floating_swap_rate;
    }
    print("avg_float is ", avg_float/double(nPaths));
    return res/double(nPaths);
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
