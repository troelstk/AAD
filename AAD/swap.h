//
//  swap.h
//  AAD
//
//  Created by Troels Tang Karlsen on 16/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef swap_h
#define swap_h

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

template<class T> T swap_price(T r_fix, T r_short, T notional,
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
