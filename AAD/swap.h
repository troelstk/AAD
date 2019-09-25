//
//  swap.h
//  AAD
//
//  Created by Troels Tang Karlsen on 16/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef swap_h
#define swap_h

#include "utilities.h"

using namespace std;

template<class T> T swap_rate(T r, T t, T Ta, T Tb,
                         vector<T> params,
                         T (*P)(T r_, T t_, T T_, vector<T> params_), T yearly_payments)
{
    T sum(0.0);
    
    T dt(1/yearly_payments);
    
    T no_payments( (Tb - Ta) * yearly_payments + 1);
    
    for(int i=1; i<no_payments; ++i){
        T Ti(Ta + double(i) * dt);
        sum += dt * P(r, t, Ti, params);
    }
    T swap_rate_res = (P(r, t, Ta, params) - P(r, t, Tb, params))/sum;
    
    return swap_rate_res;
}

template<class T> T yield_from_df(T disc, T yearly_payments, T t1, T t2){
    // Compute yield from Ta to Ti, based on discount factor. Eq. 1.13 in Brigo
    return yearly_payments * (1.0/pow(disc, 1.0/(yearly_payments * (t2 - t1))) - 1.0 );
}

template<class T> T DF_from_F(T forward, T yearly_payments){
    // Also known as FP_j(t) in Brigo
    return 1.0/(1.0 + 1.0/yearly_payments * forward);
}

template<class T> T DF_from_F(vector<T> forwards, T yearly_payments, T Ta, T Tb){
    // Computes discount factor from Ta to Tb, given forward rates from Ta to Tb
    T prod(1.0);
       
    // Loop from time alpha + 1, which is 0 here
    for(int i = 0; i<(Tb-Ta)*yearly_payments; ++i ) {
       prod *= DF_from_F(forwards[i], yearly_payments);
        //print("i is ", i, " df is ", DF_from_F(forwards[i], yearly_payments) );
    }
    
    return prod;
}



template<class T> T SR_from_F(vector<T> forwards, T yearly_payments, T Ta, T Tb){
    // Compute swap rate from Ta to Ti at time t, based on forward rates. Eq. 6.33 in Brigo
    T prod(1.0);
    T prod2(1.0);
    T sum(0.0);
    
    // Loop from time alpha + 1, which is 0 here
    for(int i = 0; i<(Tb-Ta)*yearly_payments; ++i ) {
        prod *= DF_from_F(forwards[i], yearly_payments);
    }
    
    for(int i = 0; i<(Tb-Ta)*yearly_payments; ++i ) {
        prod2 = 1.0;
        for(int j = 0; j<i+1; ++j) {
            prod2 *= DF_from_F(forwards[j], yearly_payments);
            //print(j, " ", prod2);
        }
        sum += 1.0/yearly_payments * prod2;
    }
    
    return (1.0 - prod) / sum;
}




template<class T> T vasicek_swap_price(T r_fix, T r_short, T notional,
                               T t, T Ta, T Tb, vector<T> params,
                               T (*P)(T r_, T t_, T T_, vector<T> params_), T yearly_payments)
{
    // Computes time t value of swap where term structure is given by vasicek
    T res(0.0);
   
    T no_payments( (Tb - Ta) * yearly_payments + 1);
    
    T dt( (Tb - Ta)/ (no_payments - 1) ); // Alt: 1/yearly_payments
    
    T r_float(0.0), disc(0.0);
    
    for(int i = 1; i<no_payments; ++i){
        T Ti(Ta + double(i) * dt);
        // Discounting from Ti back to Ta
        disc = P(r_short, Ta, Ti, params);
        // Compute yield from Ta to Ti
        r_float = yield_from_df(disc, yearly_payments, Ta, Ti);
        // Compute time t value of swap
        res += P(r_short, t, Ti, params) * (r_fix - r_float) * notional; // Receiver: r_fix - r_float, payer: r_float - r_fix
    }
    
    return res;
}



#endif /* swap_h */
