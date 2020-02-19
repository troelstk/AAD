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
#include "gaussians.h"

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

template<class T> T DF_from_F(T forward, double yearly_payments){
    // Also known as FP_j(t) in Brigo
    return 1.0/(1.0 + 1.0/yearly_payments * forward);
}

template<class T> T DF_from_F(vector<T> & forwards, double yearly_payments, int entry1, int entry2){
    // Computes discount factor from Ta to Tb, given forward rates from Ta to Tb
    T prod(1.0);
       
    // Loop from time alpha + 1, which is i=0 here to (beta)
    for(int i = entry1; i<entry2; ++i ) {
        prod *= DF_from_F(forwards[i], yearly_payments);
    }
    return prod;
}

template<class T> T C_ab(vector<T> & forwards, double yearly_payments, int t, int Ta, int Tb){
    // Computes sum of discount factors P(t,T_i) for T_i = Ta,..,Tb, given forward rates from t to Tb
    T sum(0.0);
    for(int i = Ta+1; i<=Tb; ++i ) {
        sum += DF_from_F(forwards, yearly_payments, t, i );
    }
    return sum;
}

template<class T> T w(vector<T> & forwards, double yearly_payments, T inp_i, int start_idx, int end_idx){
    T sum(0.0);
    
    T top_prod(1.0);
    for(int i = start_idx; i<inp_i + 1; ++i ) {
        top_prod *= DF_from_F(forwards[i], yearly_payments);
    }

    T bot_prod(1.0);
    for(int k = start_idx; k<end_idx; ++k ) {
        bot_prod = 1.0;
        for(int j = start_idx; j<k+1; ++j ) {
           bot_prod *= DF_from_F(forwards[j], yearly_payments);
        }
        sum += bot_prod;
    }
    return top_prod/sum;
}

template<class T> T vol_TFM(vector<T> & forwards, double yearly_payments, T Ta,
                            vector<vector<T>> & corr, vector<vector<T>> & vol, T swap_rate, int start_idx, int end_idx){
    // Rebonatos formula, 6.67 in Brigo. Works for time independent vol only, as specification in Table 3
    T res(0.0);
    
    for(int i = start_idx; i<end_idx; ++i){ // F10,...,F19
        for(int j = start_idx; j<end_idx; ++j){
            res += w(forwards, yearly_payments, double(i), start_idx, end_idx )
                 * w(forwards, yearly_payments, double(j), start_idx, end_idx )
            * forwards[i] * forwards[j] * corr[i][j] / swap_rate / swap_rate * Ta * vol[i][i] * vol[j][j];
        }
    }
    return res;
}

template<class T> T BlackCall(T K, T F0, T vol)
{
    // Blacks formula as defined in chapter 6.4 in Brigo
    T res(0.0);

    T d1, d2;
    
    d1 = (log(F0/K)+0.5*vol*vol)/vol;
    d2 = (log(F0/K)-0.5*vol*vol)/vol;

    res = F0 * normalCdf(d1) - K * normalCdf(d2);
    
    return res;
}


template<class T> T BlackiVol(T K, T F0, T price)
{
    // Returns the implied volatility in the Black model given K, F and a price from the Black formula
    T ub = 2;
    T lb = 0.0001;
    T mb = 0.0;
    T bs = 0.0;
    while (fabs(ub - lb) > 0.00001)
    {
        mb = 0.5*(ub + lb);
        bs = BlackCall(K, F0, mb);
        if (bs > price)
        {
            ub = mb;
        }
        else
        {
            lb = mb;
        }
    }

    return 0.5*(ub + lb);
}



template<class T> T SR_from_F(vector<T> & forwards, double yearly_payments, int entry_first, int entry_last){
    // Compute swap rate from forward rates. Eq. 6.33 in Brigo.
    T prod(1.0);
    T prod2(1.0);
    T sum(0.0);
    // Loop from time alpha + 1 to beta
    for(int i = entry_first; i < entry_last; ++i ) {
        prod *= DF_from_F(forwards[i], yearly_payments);
    }
    
    for(int i = entry_first; i<entry_last; ++i ) {
        prod2 = 1.0;
        for(int j = entry_first; j<i+1; ++j) {
            prod2 *= DF_from_F(forwards[j], yearly_payments);
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
