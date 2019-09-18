//
//  vasicek.h
//  AAD
//
//  Created by Troels Tang Karlsen on 13/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef vasicek_h
#define vasicek_h

#include <iostream>
#include "mrg32k3a.h"
#include <cmath>
#include "swap.h"


template<class T> T B(T t, T TT, T k){
    T B = 1/k * (1 - exp(-k * (TT - t)));
    return B;
}

template<class T> T A(T t, T TT, T k, T theta, T sig){
    T A = exp( (theta - sig * sig/(2 * k * k)) * (B(t, TT, k) - TT + t)
              - sig * sig/(4 * k) * B(t, TT, k) * B(t, TT, k) );
    return A;
}

// Computes discount factor from short rate, Vasicek model
template<class T> T P(T r, T t, T TT, vector<T> params){
    T k(params[0]), theta(params[1]), sig(params[2]);
    T P = A(t, TT, k, theta, sig) * exp( -B(t, TT, k) * r ) ;
    return P;
}

template<class T> T vasicek_swap(vector<T> params, T t, T Ta, T Tb, T r_fix, T notional, T r0,
                            int seed1, int seed2, int nPaths, int nSteps, T yearly_payments)
{
    // Simulates rates from Vasicek model from time t to time Ta and
    //  computes swap price from Ta to Tb, with 
    T k(params[0]), theta(params[1]), sig(params[2]);
    T r_sim(r0);
    
    T dt((Tb-Ta)/double(nSteps));
    
    T res = 0.0;
    
    mrg32k3a myRNG(seed1, seed2);
    myRNG.init(nSteps);
    vector<double> gauss(nSteps);
    
    for(int i=0; i<nPaths; ++i){
        myRNG.nextG(gauss);
        r_sim = r0;
        for(int j=0; j<nSteps; ++j){
            r_sim += k * (theta - r_sim) * dt + sig * sqrt(dt) * gauss[j];
        }
        //cout << "sim rate is " << r_sim << endl;
        T temp = swap_price(r_fix, r_sim, notional, t, Ta, Tb, params, P, yearly_payments);
        res += temp;
        //cout << temp << endl;
    }
    
    res /= nPaths;
    
    return res;
}


#endif /* vasicek_h */
