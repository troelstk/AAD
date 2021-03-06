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
//#include "Heston.h"
#include "lsm.h"
//#include "vasicek.h"
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
    int seed1, nSteps, nPaths, nPaths_presim, M, dim_n, seed2, nSteps_y;

    seed1 = 59;
    nSteps = 4;
    nPaths = 100;
    
    t = 0.0;
    Ta = 1.0; // Payment at time 2, 3 and 4
    Tb = 4.0;
    int int_Ta = int(Ta), int_Tb = int(Tb), int_t = int(t);
    /*
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
    
    double Cab = C_ab( F, yearly_payments, int_t, int_Ta, int_Tb);
    //print("Cab is ", Cab);
    // Black(T K, T F0, T vol)
    double black = Cab * BlackCall( r_fix, swap_rate, rebonato);
    
    double disc = DF_from_F(F[0], Ta);
    
    print("Black price ", disc * black * notional);
    
    vector<double> exTimes = {0.0, 1.0, 2.0, 3.0};  // , 2.0, 3.0
    
    
    
    clock_t begin_timeBswap = clock();*/
    
   

    
    
    
    
    
    // Test case as in 8.2 in Brigo
    
    print("\nTEST CASE: ");
    t = 0.0;
    Ta = 9; // Start of swap
    Tb = 20; // End of swap
    
    //vector<double> exTimes20   = {0.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0};
    //vector<double> exTimes20   = {0.0, 9.0};
    
    int first_ex = Ta; // First time swap can be exercised
    int last_ex = Tb-1; // Last time swap can be exercised

    vector<double> exTimes20(last_ex-first_ex+2, 0.0);
    for(int i = 1; i< last_ex-first_ex+2; ++i){
        exTimes20[i] = double(i+first_ex-1);
    }
    print(exTimes20);
    
    
    int_Ta = int(Ta); int_Tb = int(Tb); int_t = int(t);
    
    dim_n = 2;
    yearly_payments = 1.0;
    notional = 100.0;
    vector<double> F20 = {0.0469, 0.0501, 0.0560, 0.0584, 0.0600, 0.0613, 0.0628, 0.0627, 0.0629, 0.0623, 0.0630, 0.0636, 0.0643, 0.0648, 0.0653, 0.0640, 0.0630, 0.0618, 0.0607, 0.0594 };
    /*vector<double> F20(int_Tb);
    for(int i = 0; i<Tb; ++i){
        F20[i] = 1.0/70.0*pow(i,0.05)+0.0469;
        //F20[i] = 0.10;
    }
    print(F20);*/
    
    
    M = int_Tb;
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
            corrA[i][j] = cos( Theta[i-1] - Theta[j-1]) ;
            //corrA[i][j] = 1.0;
            //corrA[i][j] = 1.0/(double((abs(j-i)*abs(i-j)+8.0*Tb))/(8*Tb));
        }
        vol20[i][i] = Phi[i] ;
        //vol20[i][i] = 0.2 ;
        //vol20[i][i] = 0.045*pow(i,0.2)*exp(-i*0.025) + 0.095 ;
    }
    //print(corrA);
    
    double rebonato2 = sqrt(vol_TFM(F20, yearly_payments, Ta, corrA, vol20, swap_rate20, int(Ta) + 1, int_Tb) );
    print("Rebonato vol is ", rebonato2/sqrt(Ta));
    double C = C_ab( F20, yearly_payments, int_t, int_Ta, int_Tb);
    print("C is ", C);
    double black20 = BlackCall( swap_rate20, swap_rate20, rebonato2);
    print("Black price ", C * notional * black20 );
    
    
    nSteps = 36;
    nPaths = 10000;
    seed1 = 3;
    
    double lmm20 = -1;
    clock_t begin_time2 = clock();
    lmm20 = LMM_swaption(vol20, corrA, F20, t, Ta, Tb, r_fix, notional, seed1, nPaths, nSteps, yearly_payments, dim_n);
    auto time2 =  float( clock () - begin_time2 )/ CLOCKS_PER_SEC;
    print("LMM EUR Swaption price is ", lmm20);
    print("Calculation time: ", time2, " seconds");
    
    double BlackImpVol = BlackiVol(swap_rate20, swap_rate20, lmm20 / ( C * notional) );
    print("Black implied vol of simulation is ", BlackImpVol/sqrt(Ta) );
    
    
    print("\nBermuda test: ");
    nPaths_presim = 2000;
    nPaths = 10000; // Main
    
    nSteps_y = 4;
    seed1 = 245;
    seed2 = 295;
    print("Seed1: ", seed1, ", seed2: ", seed2, ", Main paths: ", nPaths, " yearly steps: ", nSteps_y);
    

    
    print("\nAAD test:");
    vector<number> F_Risk, thetaRisk, PhiRisk;
    
    for(auto & x : F20)   {F_Risk.push_back(number(x));}
    
    for(auto & x : Theta) {thetaRisk.push_back(number(x));}
    for(auto & x : Phi)   {PhiRisk.push_back(number(x));}
    
    vector<vector<number>> corrRisk(M, vector<number>(M, number(0.0)));
    vector<vector<number>> volRisk(M,  vector<number>(M, number(0.0)));
    
    // Fill Vol and corr
    for(int i = 1; i<M; i++){
        for(int j = 1; j<i+1; j++){
            corrRisk[i][j] = cos( thetaRisk[i-1] - thetaRisk[j-1]);
            //corrRisk[i][j] = number(1.0/(double((abs(j-i)*abs(i-j)+8.0*Tb))/(8*Tb)));
            //corrRisk[i][j] = 1.0;
            //print("row ", i, " col ", j, " ", corrRisk[i][j].value());
        }
        //volRisk[i][i] = number(0.045*pow(i,0.2)*exp(-i*0.025) + 0.095 );
        volRisk[i][i] = PhiRisk[i];
        //volRisk[i][i] = number(0.2);
    }
    //print("other side");
    for(int i = 1; i<M; i++){
        for(int j = 1; j<i; j++){
            corrRisk[j][i] = corrRisk[i][j];
            //print("row ", j, " col ", i, " ", corrRisk[i][j].value());
        }
    }
    number fixedRateRisk(r_fix);
    
    
    
    
    print("\nAAD test Bermuda:");
    clock_t bermAADstart = clock();

    number bermudanAAD;
    //bermudanAAD = LMM_BermudaSwaptionAAD(volRisk, corrRisk, F_Risk, exTimes20, t, Ta, Tb, fixedRateRisk, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
    bermudanAAD = LMM_BermudaSwaptionAADChol(volRisk, corrRisk, F_Risk, exTimes20, t, Ta, Tb, fixedRateRisk, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
    auto timeAADBermuda = float( clock() - bermAADstart )/ float( CLOCKS_PER_SEC);
    print("Bermudan value AAD: ", bermudanAAD.value(), " in ", timeAADBermuda, " seconds");
    
    //for(auto & x : volRisk[1]){print("Vol risk is ", x.adjoint()/double(nPaths)) ;}
    //double FR_AAD = fixedRateRisk.adjoint()/double(nPaths);
    
    double eps = 10e-12;
    
    /* Print all F adjoints */
    print("F adjoints");
    for(int i=0; i<Tb; ++i) { // i<Tb
       double AAD_approx = F_Risk[i].adjoint()/double(nPaths);
       printf("%1.2f & ", AAD_approx);
    }
    cout << "\n";
    /* Print all vol adjoints */
    print("Vol adjoints");
    for(int i=0; i<Tb; ++i) { // i<Tb
       //double AAD_approx = PhiRisk[i].adjoint()/double(nPaths);
       double AAD_approx = volRisk[i][i].adjoint()/double(nPaths);
       printf("%1.2f & ", AAD_approx);
    }
    cout << "\n";
       
    /* Print all correlation adjoints */
    /*print("Correlation Adjoints");
    for(int i=9; i<Tb; ++i) { // i<Tb
        printf("%1.0d & ", i);
        for(int j=9; j<Tb; ++j) { // j<i+1
           double AAD_approx = corrRisk[i][j].adjoint()/double(nPaths);
           if( i != j ) {
               printf("%5.3f & ", AAD_approx);
           }
           else {
               printf("- & ");
           }
       }
       cout << "\n";
    }*/
    
    //number::tape->clear();
    
    // Compute unbumped value:
    clock_t begin_timeBswap2 = clock();

    double bermudan_unbumped = LMM_BermudaSwaption2(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
    auto timeBSwap3 =  float( clock() - begin_timeBswap2 )/ CLOCKS_PER_SEC;
    print("Bermudan value non-AAD: ", bermudan_unbumped, " in ", timeBSwap3, " seconds");

    print("AAD: ", float(timeAADBermuda) / float(timeBSwap3), " times slower");
    
    /* Compare forward rate adjoints with bump-and-revalue */
    /*print("Compare forward rate adjoints with bump-and-revalue");
    for(int idx=0; idx<Tb; ++idx) {
        // Bump one at the time and compute FD approx:
        F20[idx] += eps;
        double bermudan_bumped = LMM_BermudaSwaption2(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
        F20[idx] -= eps;
        double FD_approx = (bermudan_bumped - bermudan_unbumped)/eps;
        
        double AAD_approx = F_Risk[idx].adjoint()/double(nPaths);
        //printf("%11.10f & ", (AAD_approx-FD_approx )/ AAD_approx*100 );
        //printf("%1.2f %% & ", (AAD_approx-FD_approx )/ AAD_approx*100 );
        printf("%2.1d: %11.10f,%11.10f,%11.10f\n", idx, AAD_approx, FD_approx, AAD_approx/FD_approx);
        //print(float( clock() - begin_timeBswap2 )/ CLOCKS_PER_SEC);
    }
    cout << "\n";*/
    
    /* Compare vol adjoints with bump-and-revalue */
    /*print("Compare vol adjoints with bump-and-revalue");
    for(int idx=9; idx<Tb; ++idx) {
        // Bump one at the time and compute FD approx:
        vol20[idx][idx] += eps;
        double bermudan_bumped = LMM_BermudaSwaption2(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
        vol20[idx][idx] -= eps;
        double FD_approx = (bermudan_bumped - bermudan_unbumped)/eps;
        
        double AAD_approx = PhiRisk[idx].adjoint()/double(nPaths);
        //printf("%1.2f %% & ", (AAD_approx-FD_approx )/ AAD_approx*100 );
        printf("%2.1d: %11.10f,%11.10f,%11.10f\n", idx, AAD_approx, FD_approx, AAD_approx/FD_approx);
        //print(float( clock() - begin_timeBswap2 )/ CLOCKS_PER_SEC);
    }
    cout << "\n";*/
    
    /* Compare correlation adjoints with bump-and-revalue */
    print("Compare correlation adjoints with bump and revalue");
    for(int i=9; i<int_Tb; ++i) { // i<Tb
        for(int j=9; j<int_Tb; ++j) { // j<i+1
        // Bump one at the time and compute FD approx:
            if(i != j) {
                corrA[i][j] += eps;
                corrA[j][i] += eps;
                double bermudan_bumped = LMM_BermudaSwaption2(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix, notional, seed1, seed2, nPaths, nPaths_presim, nSteps_y, yearly_payments, dim_n);
                corrA[i][j] -= eps;
                corrA[j][i] -= eps;
                
                double FD_approx = (bermudan_bumped - bermudan_unbumped)/eps;
                
                double AAD_approx = corrRisk[i][j].adjoint()/double(nPaths);
                printf("%2.1d, %2.1d: %11.10f,%11.10f,%11.10f\n", i, j, AAD_approx, FD_approx, AAD_approx/FD_approx);
                //printf("%1.4f & ", AAD_approx/ FD_approx );
                //print(float( clock() - begin_timeBswap2 )/ CLOCKS_PER_SEC);
            }
            else {
                //printf("%1.4f, ", - );
            }
        }
        cout << "\n";
    }
    
    print("Bump and revalue: ", float( clock() - begin_timeBswap2 )/ CLOCKS_PER_SEC, " seconds");
    
   
    

    
    /*print("\nAAD test European:");
    nPaths = 10000;
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
    
    /* Test number of paths needed in pre sim:
     vector<int> nPath_vec = {100, 200, 400, 800, 1000, 1500, 2000, 3000, 4000, 5000, 10000, 20000};
    
    for(auto & x : nPath_vec){
        clock_t time1 = clock();
        double bermudan_paths = LMM_BermudaSwaption(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix, notional, seed1, seed2, 70000, x, nSteps_y, yearly_payments, dim_n);
        printf("%11.10f,%11.10f\n", bermudan_paths, float( clock () - time1 )/ CLOCKS_PER_SEC );
    }*/
    
    /* Test number of paths needed in main sim:
     vector<int> nPath_vec = {30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000};
    for(auto & x : nPath_vec){
        clock_t time1 = clock();
        double bermudan_paths = LMM_BermudaSwaption(vol20, corrA, F20, exTimes20, t, Ta, Tb, r_fix, notional, seed1, seed2, x, 3000, nSteps_y, yearly_payments, dim_n);
        printf("%11.10f,%11.10f\n", bermudan_paths, float( clock () - time1 )/ CLOCKS_PER_SEC );
    }*/
    
    return 0;
}

