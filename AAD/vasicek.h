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
#include <cmath>
#include <cstdlib>

#include "mrg32k3a.h"
#include "swap.h"
#include "gaussians.h"
#include "AAD.h"


template<class T> T B(T t, T TT, T k){
    // B function in Vasicek model
    T B = 1.0/k * (1.0 - exp(-k * (TT - t)));
    return B;
}

template<class T> T A(T t, T TT, T k, T theta, T sig){
    // A function in Vasicek model
    T A = exp( (theta - sig * sig/(2.0 * k * k)) * (B(t, TT, k) - TT + t)
              - sig * sig/(4.0 * k) * B(t, TT, k) * B(t, TT, k) );
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
    //  computes price of swap starting at Ta and lasting until Tb
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
        T temp = vasicek_swap_price(r_fix, r_sim, notional, t, Ta, Tb, params, P, yearly_payments);
        res += temp;
    }
    
    res /= nPaths;
    
    return res;
}

template<class T> T vasicek_swaption(vector<T> params, T t, T Ta, T Tb, T r_fix, T notional, T r0,
                                 int seed1, int seed2, int nPaths, int nSteps, T yearly_payments, T strike)
{
    // Computes Vasicek european swaption price with maturity Ta and tenor Tb
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
        // Compute time Ta price of a swap from Ta+1 to Tb
        T price = vasicek_swap_price(r_fix, r_sim, notional, Ta, Ta, Tb, params, P, yearly_payments);
        res += max(price - strike, 0.0);
    }
    // Discount result back from time Ta to time t
    res = P(r_sim, t, Ta, params) * res / nPaths;
    
    
    return res;
}

template<class T> T vasicek_swaption_aad(vector<T> params, T t, T Ta, T Tb, T r_fix, T notional, T r0,
                                     int seed1, int seed2, int nPaths, int nSteps, T yearly_payments, T strike)
{
    // Computes Vasicek european swaption price with maturity Ta and tenor Tb
    T k(params[0]), theta(params[1]), sig(params[2]);
    T r_sim(r0);
    
    T dt((Tb-Ta)/double(nSteps));

    T res(0.0);
    
    mrg32k3a myRNG(seed1, seed2);
    myRNG.init(nSteps);
    vector<double> gauss(nSteps);
    
    number::tape->mark();
    
    for(int i=0; i<nPaths; ++i){
        number::tape->rewindToMark();
        myRNG.nextG(gauss);
        r_sim = r0;
        for(int j=0; j<nSteps; ++j){
            r_sim += k * (theta - r_sim) * dt + sig * sqrt(dt) * gauss[j];
        }
        // Compute time Ta price of a swap from Ta+1 to Tb
        T price = vasicek_swap_price(r_fix, r_sim, notional, Ta, Ta, Tb, params, P, yearly_payments);
        
        T payoff(max(price - strike, 0.0));
        
        payoff.propagateToMark();
        
        res += payoff;
    }
    // Discount result back from time Ta to time t
    res = P(r_sim, t, Ta, params) * res / double(nPaths);
    number::propagateMarkToStart();

    return res;
}

void test_vasicek()
{
    double mat = 4.0;
    int seed1 = 12, seed2 = 2344, nPaths = 10000;
    int nSteps = 25;
    double r0 = 0.025, t = 0.0, yearly = 1.0 ;
    
    double k = 0.2; // Speed of mean reversion
    double theta = 0.04; // Long term average rate
    double sig = 0.02; // instantaneous volatility
    
    double r_fix = 0.04, r_short = 0.02, notional = 100, Ta = 1.0, Tb = 4.0, yearly_payments = 1.0;
    double strike = 0.0;
    
    vector<double> params = {k, theta, sig};
    
    double swap_rate_1 = swap_rate(r0, t, t, mat, params, P, yearly);
    cout << "Swap rate is " << swap_rate_1 << endl;
    
    
    double swap_price1 = vasicek_swap_price(r_fix, r_short, notional, t, Ta, Tb, params, P, yearly_payments);
    
    cout << "Swap price is " << swap_price1 << endl;
    
    double sim_swap = vasicek_swap(params, t, Ta, Tb, r_fix, notional, r0,
                                   seed1, seed2, nPaths, nSteps, yearly_payments);
    cout << "simulated swap value is " << sim_swap << endl;
    

    nPaths = 10000;
    double sim_swaption = vasicek_swaption(params, t, Ta, Tb, r_fix, notional, r0,
                                           seed1, seed2, nPaths, nSteps, yearly_payments, strike);
    cout << "simulated swaption value is " << sim_swaption << endl;
    
    double bump = 0.00001;
    params[1] += bump;
    double sim_swaption_bump = vasicek_swaption(params, t, Ta, Tb, r_fix, notional, r0,
                                                seed1, seed2, nPaths, nSteps, yearly_payments, strike);
    cout << "FD theta derivative is " << (sim_swaption_bump-sim_swaption)/bump << endl;
    
    params[1] -= bump;
    params[0] += bump;
    sim_swaption_bump = vasicek_swaption(params, t, Ta, Tb, r_fix, notional, r0,
                                         seed1, seed2, nPaths, nSteps, yearly_payments, strike);
    cout << "FD k derivative is " << (sim_swaption_bump-sim_swaption)/bump << endl;
    
    
    
    number::tape->rewind();
    number k1(k); // Speed of mean reversion
    number theta1(theta); // Long term average rate
    number sig1(sig); // instantaneous volatility
    vector<number> params1 = {k1, theta1, sig1};
    number t1(t), Ta1(Ta), Tb1(Tb), r_fix1(r_fix), notional1(notional), r01(r0);
    number yearly_payments1(yearly_payments), strike1(strike);
    
    number sim_swaption_aad = vasicek_swaption_aad(params1, t1, Ta1, Tb1, r_fix1, notional1, r01,
                                                   seed1, seed2, nPaths, nSteps, yearly_payments1, strike1);
    cout << "simulated swaption value with aad is " << sim_swaption_aad.value() << endl;
    number::tape->rewind();
    cout << "Theta deriv is " << theta1.adjoint()/double(nPaths) << endl;
    cout << "k deriv is " << k1.adjoint()/double(nPaths) << endl;

}


#endif /* vasicek_h */
