//
//  main.cpp
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//


//#include "/usr/local/Cellar/armadillo/9.700.2/include/armadillo"


//using namespace arma;

#include <iostream>
#include <vector>
#include <ctime>

#include "AAD.h"
#include "Heston.h"
#include "lsm.h"
#include "vasicek.h"
#include "swap.h"
#include "LMM.h"
#include "utilities.h"
#include "Bermuda.h"
#include "BermudaAAD.h"

using namespace std;



int main(int argc, const char * argv[]) {
   
    //test_heston();
    //test_vasicek();
    double t , notional, Ta, Tb, yearly_payments,  r_fix;
    int seed1, nSteps, nPaths, M, dim_n, seed2;
    
    
    seed1 = 59;
    dim_n = 2;
    nSteps = 4;
    nPaths = 100;
    
    t = 0.0;
    Ta = 1.0; // Payment at time 2, 3 and 4
    Tb = 4.0;
    int int_Ta = int(Ta), int_Tb = int(Tb), int_t = int(t);
    yearly_payments = 1.0;
    notional = 100.0;
    
    M = 4; // Dimension of forward curve
    vector<vector<double>> vol(M,  vector<double>(M));
    vector<vector<double>> corr(M, vector<double>(M));
    vector<double> F = {0.01, 0.02, 0.03, 0.04}; // F_1, F_2, F_3, F_4
    double swap_rate = SR_from_F(F, yearly_payments, int_Ta, int_Tb);
    print( "Swap rate is ", swap_rate );
    r_fix = swap_rate;
    //vector<double> r_fix_vec = {0.05, 0.01, 0.02, swap_rate, 0.035, 0.04, 0.045};
    
    // Table 3: Constant vol for each F regardless of time t
    vol[0] = {0.20, 0.00, 0.00, 0.00};
    vol[1] = {0.25, 0.25, 0.00, 0.00};
    vol[2] = {0.30, 0.30, 0.30, 0.00};
    vol[3] = {0.30, 0.30, 0.30, 0.30};
    // Burde bruge table 5 med time dependent vol
    
    corr[0] = {1.0, 0.8, 0.7, 0.6};
    corr[1] = {0.8, 1.0, 0.8, 0.7};
    corr[2] = {0.7, 0.8, 1.0, 0.8};
    corr[3] = {0.6, 0.7, 0.8, 1.0};
    // Bør lave mere sofistikeret corr
    
    clock_t begin_time = clock();
    
    double lmm = LMM_swaption(vol, corr, F, t, Ta, Tb, r_fix, notional, seed1, nPaths, nSteps, yearly_payments, dim_n);
    auto time =  float( clock () - begin_time )/ CLOCKS_PER_SEC;
    print("LMM Swaption price is ", lmm);
    print("Calculation time: ", time, " seconds");

    
    double sum = 0.0;
    for(int i = int(Ta); i<int(Tb); ++i){
        double weight = w(F, yearly_payments, double(i), int(Ta));
        //print("weight is ", weight);
        sum += weight * F[i];
    }
    print("swap rate from weights is ", sum);
    
    double rebonato = sqrt(vol_TFM(F, yearly_payments, Ta, corr, vol, swap_rate, Ta) );
    print("Rebonato vol is ", rebonato);
    
    double Cab = C_ab( F, yearly_payments, int_Ta, int_Tb);
    //print("Cab is ", Cab);
    // Black(T K, T F0, T vol)
    double black = Cab * BlackCall( r_fix, swap_rate, rebonato);
    
    double disc = DF_from_F(F[0], Ta);
    
    print("Black price ", disc * black * notional);
    
    vector<double> exTimes = {0.0, 1.0, 2.0, 3.0};  // , 2.0, 3.0
    
    int nPaths_presim = 1000;
    nPaths = 100000;
    int nSteps_y = 4;
    
    clock_t begin_timeBswap = clock();
    
    //double bermudan = LMM_BermudaSwaption(vol, corr, F, exTimes, t, Ta, Tb, r_fix, notional, seed1, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
    //auto timeBSwap =  float( clock () - begin_timeBswap )/ CLOCKS_PER_SEC;
    
    //print("Bermuda price is ", bermudan, ", time ", timeBSwap, " seconds");
    

    
    
    
    
    
    // Test case as in 8.2 in Brigo
    
    print("\nTEST CASE: ");
    // T0 = 1y = 1.0
    t = 0.0;
    Ta = 9;
    Tb = 20;
    int_Ta = int(Ta); int_Tb = int(Tb); int_t = int(t);
    
    dim_n = 2;
    yearly_payments = 1.0;
    notional = 100.0;
    vector<double> F20 = {0.0469, 0.0501, 0.0560, 0.0584, 0.0600, 0.0613, 0.0628, 0.0627, 0.0629, 0.0623,
                          0.0630, 0.0636, 0.0643, 0.0648, 0.0653, 0.0640, 0.0630, 0.0618, 0.0607, 0.0594 };
    
    // Test 1.b, formulation 7 with a=0, b=0, c=1, d=0, all correlations set to 1
    M = 20;
    vector<vector<double>> vol20(M, vector<double>(M));
    vector<vector<double>> corrA(M, vector<double>(M));

    double swap_rate20 = SR_from_F(F20, yearly_payments, int_Ta, int_Tb); // 9'th entry is F10,
    print( "Swap rate is ", swap_rate20 );
    r_fix = swap_rate20;
    
    vector<double> Phi = {
        -1000, 0.149, 0.159, 0.153, 0.1450, 0.1360, 0.1270, 0.1210, 0.1180, 0.1140,
        0.111, 0.108, 0.105, 0.102, 0.0989, 0.0978, 0.0974, 0.0969, 0.0965, 0.0961};
    vector<double> Theta = {
        0.0147, 0.0643, 0.1032, 0.1502, 0.1969, 0.2239, 0.2771, 0.2950,
        0.3630, 0.3810, 0.4217, 0.4836, 0.5204, 0.5418, 0.5791, 0.6496,
        0.6679, 0.7126, 0.7659};

    // Fill Vol and corr
    for(int i = 1; i<M; i++){
        for(int j = 1; j<M; j++){
            vol20[i][i] = Phi[i] ;
            corrA[i][j] = cos( Theta[i-1] - Theta[j-1]) ;
        }
    }
    //print(corrA);

    double rebonato2 = sqrt(vol_TFM(F20, yearly_payments, Ta, corrA, vol20, swap_rate20, int(Ta) + 1) );
    print("Rebonato vol is ", rebonato2/sqrt(Ta));
    
    double C = C_ab( F20, yearly_payments, int_Ta, int_Tb);
    print("C is ", C);
    // Black(T K, T F0, T vol)
    
    double disc2 = DF_from_F(F20, yearly_payments, int_t, int_Ta);
    
    double black20 = BlackCall( swap_rate20, swap_rate20, rebonato2);
    print("Black price ", disc2 * C * notional * black20 );
    double BlackImpVol3 = BlackiVol(swap_rate20, swap_rate20, black20 );
    print("Black20 implied vol is ", BlackImpVol3/sqrt(Ta) );
    
    
    nSteps = 36;
    nPaths = 10000;
    seed1 = 3;
    
    double lmm20 = -1.0;
    clock_t begin_time2 = clock();
    lmm20 = LMM_swaption(vol20, corrA, F20, t, Ta, Tb, r_fix, notional, seed1, nPaths, nSteps, yearly_payments, dim_n);
    auto time2 =  float( clock () - begin_time2 )/ CLOCKS_PER_SEC;
    print("LMM EUR Swaption price is ", lmm20);
    print("Calculation time: ", time2, " seconds");
    
    double BlackImpVol = BlackiVol(swap_rate20, swap_rate20, lmm20 / (disc2 * C * notional) );
    print("Black implied vol of simulation is ", BlackImpVol/sqrt(Ta) );
    
    
    print("\nBermuda test: ");
    nPaths_presim = 1000;
    nPaths = 1000; // Main
    
    nSteps_y = 4;
    seed2 = 41232;
    seed1 = 3;
    
    clock_t begin_timeBswap2 = clock();
    vector<double> exTimes20 = {0.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0};
    //vector<double> exTimesEur = {0.0, 9.0};
    int vol_idx = 15;
    //print("vol ", vol20[vol_idx][vol_idx] );
    
    double eps = 0.000001;
    double bermudan2 = LMM_BermudaSwaption(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
    print("Bermudan price ", bermudan2);
    vol20[vol_idx][vol_idx] += eps;
    //print("vol ", vol20[vol_idx][vol_idx] );
    double bermudan3 = LMM_BermudaSwaption(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
    vol20[vol_idx][vol_idx] -= eps;
    double bermudan4 = LMM_BermudaSwaption(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix+eps, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
    
    auto timeBSwap3 =  float( clock () - begin_timeBswap2 )/ CLOCKS_PER_SEC / 3.0;

    double vol_FD = (bermudan3 - bermudan2)/eps;
    double FR_FD = (bermudan4 - bermudan2)/eps;
    print("Vol FD: ", vol_FD);
    print("Fixed rate FD: ", FR_FD, ", time is ", timeBSwap3);
    
    
    print("\nAAD test:");
    
    vector<number> F_Risk, thetaRisk, PhiRisk;
    
    //for(auto & x : exTimes20) {exTimesRisk.push_back(number(x));}
    for(auto & x : F20) {F_Risk.push_back(number(x));}
    for(auto & x : Theta) {thetaRisk.push_back(number(x));}
    for(auto & x : Phi) {PhiRisk.push_back(number(x));}
    
    vector<vector<number>> corrRisk(M, vector<number>(M));
    vector<vector<number>> volRisk(M, vector<number>(M));
    
    // Fill Vol and corr
    for(int i = 1; i<M; i++){
        for(int j = 1; j<M; j++){
            volRisk[i][j] =  PhiRisk[i];
            corrRisk[i][j] = cos( thetaRisk[i-1] - thetaRisk[j-1]);
        }
    }
    
    number fixedRateRisk(r_fix);
    double yearly = 1.0;
    
    
    /*
    print("\nAAD test European:");
    nPaths = 1000;
    clock_t eurAADstart = clock();
    number EUR_AAD = LMM_swaptionAAD(volRisk, corrRisk, F_Risk,
                                     t, Ta, Tb, fixedRateRisk, notional,
                                     seed1, nPaths, nSteps, yearly_payments, dim_n) ;
    print( "EUR AAD took ",  float( clock() - eurAADstart )/ CLOCKS_PER_SEC);
    print("European value AAD: ", EUR_AAD.value());
    print("Vol adjoint is ", volRisk[vol_idx][vol_idx].adjoint()/double(nPaths),
          " val is ", volRisk[vol_idx][vol_idx].value());
    print("Fixed rate adjoint is ", fixedRateRisk.adjoint()/double(nPaths) );
    
    
    eps = 0.000001;
    double lmm1, lmm2;
    vol20[vol_idx][vol_idx] += eps;
    lmm1 = LMM_swaption(vol20, corrA, F20, t, Ta, Tb, r_fix,
                        notional, seed1, nPaths, nSteps, yearly_payments, dim_n);
    vol20[vol_idx][vol_idx] -= eps;

    lmm2 = LMM_swaption(vol20, corrA, F20, t, Ta, Tb, r_fix,
                        notional, seed1, nPaths, nSteps, yearly_payments, dim_n);
    print("FD approx vol", vol_idx, " is ", (lmm1 - lmm2)/eps);
    
    print("Ratio vol ", 1.0 /( volRisk[vol_idx][vol_idx].adjoint()/double(nPaths) / ((lmm1 - lmm2)/eps) ) );

    lmm1 = LMM_swaption(vol20, corrA, F20, t, Ta, Tb, r_fix + eps,
                        notional, seed1, nPaths, nSteps, yearly_payments, dim_n);
    //print("fixed ", r_fix );
    print("FD approx fixed rate: ", (lmm1 - lmm2)/eps);
    
    print("Ratio r_fix ", fixedRateRisk.adjoint()/double(nPaths) / ((lmm1 - lmm2)/eps) );
    
    number::tape->clear();
    
    double blackBump = BlackCall( swap_rate20 + eps, swap_rate20, rebonato2);
    print("Black bumped price ", disc2 * C * notional * blackBump );
    print("Black r_fix deriv: ", disc2 * C * notional * (blackBump - black20)/eps);*/
    
    print("\nAAD test Bermuda:");
    clock_t bermAADstart = clock();
    number bermudanAAD;
    bermudanAAD = LMM_BermudaSwaptionAAD(volRisk, corrRisk, F_Risk, exTimes20, t, Ta, Tb, fixedRateRisk, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly, dim_n);
    auto timeAADBermuda = float( clock() - bermAADstart )/ float( CLOCKS_PER_SEC);
    print("Bermudan value AAD: ", bermudanAAD.value(), " in ", timeAADBermuda, " seconds");
    
    //for(auto & x : volRisk[1]){print("Vol risk is ", x.adjoint()/double(nPaths)) ;}
    double FR_AAD = fixedRateRisk.adjoint()/double(nPaths);
    double vol_AAD = volRisk[vol_idx][vol_idx].adjoint()/double(nPaths);
    print("Vol adjoint is ", vol_AAD);
    print("Fixed rate adjoint is ", FR_AAD);
    
    print("Vol ratio ", vol_FD / vol_AAD);
    print("Fixed rate ratio ", FR_FD / FR_AAD);
    
    print("\nAAD: ", float(timeAADBermuda) / float(timeBSwap3), " times slower");
    
    return 0;
}


