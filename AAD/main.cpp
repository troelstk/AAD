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

using namespace std;



int main(int argc, const char * argv[]) {
   
    //test_heston();
    //test_vasicek();
    
    /*clock_t begin_time = clock();
    lsm(200000, 50);
    auto time_aad =  float( clock () - begin_time )/  CLOCKS_PER_SEC;
    cout << "Calculation time: " << time_aad  << " seconds" << endl;
    double res = bsMonteCarlo(100, 0.2, 1, 100, 123, 234);
    cout << "Result is " <<  res <<  endl;*/

    
    int seed1 = 434;
    int nSteps = 4, nPaths = 100000;
    double t = 0.0;
    
    double Ta = 1.0, Tb = 4.0, yearly_payments = 1.0;
    double notional = 100.0, r_fix;
    
    int M = 4; // Dimension of forward curve
    vector<vector<double>> vol(M,  vector<double>(M));
    vector<vector<double>> corr(M, vector<double>(M));
    vector<double> F = {0.01, 0.02, 0.03, 0.04};
    double swap_rate = SR_from_F(F, yearly_payments, int(Ta), int(Tb));
    print( "Swap rate is ", swap_rate );
    r_fix = swap_rate;
    vector<double> r_fix_vec = {0.05, 0.01, 0.02, swap_rate, 0.035, 0.04, 0.045};
    
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

    double lmm = LMM_swaption(vol, corr, F, t, Ta, Tb, r_fix, notional, seed1, nPaths, nSteps, yearly_payments);
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
    
    double Cab = C_ab( F, yearly_payments, int(Ta), int(Tb));
    // Black(T K, T F0, T vol)
    double black = Cab * BlackCall( r_fix, swap_rate, rebonato);
    
    double disc = DF_from_F(F[0], 1.0);
    
    //print("Disc is ", disc);
    print("Black price ", disc * black * notional);

    
    
    
    
    
    
    // Test case as in 8.2 in Brigo
    
    print("\nTEST CASE: ");
    Ta = 9;
    Tb = 20;
    vector<double> F20 = {0.0469, 0.0501, 0.0560, 0.0584, 0.0600, 0.0613, 0.0628, 0.0627, 0.0629, 0.0623,
                          0.0630, 0.0636, 0.0643, 0.0648, 0.0653, 0.0640, 0.0630, 0.0618, 0.0607, 0.0594 };
    
    // Test 1.b, formulation 7 with a=0, b=0, c=1, d=0, all correlations set to 1
    M = 20;
    vector<vector<double>> vol20(M,  vector<double>(M));
    vector<vector<double>> corr20(M, vector<double>(M));

    double swap_rate20 = SR_from_F(F20, yearly_payments, int(Ta)+1, int(Tb));
    print( "Swap rate is ", swap_rate20 );
    r_fix = swap_rate20;
    
    vector<double> Phi = {
        0.000, 0.149, 0.159, 0.153, 0.1450, 0.1360, 0.1270, 0.1210, 0.1180, 0.1140,
        0.111, 0.108, 0.105, 0.102, 0.0989, 0.0978, 0.0974, 0.0969, 0.0965, 0.0961};
    // Fill Vol and corr
    for(int i = 0; i<M; i++){
        for(int j = 0; j<M; j++){
            corr20[i][j] = 1.0;
            vol20[i][j] = (i == j) & (i>0) & (j>0) ? Phi[i] : 0.0;
        }
    }
    double sum2 = 0.0;
    for(int i = int(Ta)+1; i<int(Tb); ++i){
        double weight = w(F20, yearly_payments, double(i), int(Ta)+1);
        sum2 += weight * F20[i];
    }
    print("swap rate from weights is ", sum2);

    double rebonato2 = sqrt(vol_TFM(F20, yearly_payments, Ta, corr20, vol20, swap_rate20, int(Ta) + 1) );
    print("Rebonato vol is ", rebonato2/sqrt(Ta));
    
    double C = C_ab( F20, yearly_payments, int(Ta)+1, int(Tb));
    
    // Black(T K, T F0, T vol)
    disc = DF_from_F(F20, yearly_payments, t, Ta);
    print("Disc is ", disc);
    
    double black20 = BlackCall( swap_rate20, swap_rate20, rebonato2);
    print("Black price ", disc * C * notional * black20 );
    
    //mat X = randn<mat>(5,5);
    //mat Y = X.t()*X;
    //mat R1 = chol(Y);
    //mat R2 = chol(Y, "lower");

    nSteps = 20;
    nPaths = 20000;
    
    double lmm20 = -1.0;


    lmm20 = LMM_swaption(vol20, corr20, F20, t, Ta, Tb, r_fix, notional,
                          seed1, nPaths, nSteps, yearly_payments);
    time =  float( clock () - begin_time )/ CLOCKS_PER_SEC;
    print("LMM Swaption price is ", lmm20);
    print("Calculation time: ", time, " seconds");
    
    double BlackImpVol = BlackiVol(swap_rate20, swap_rate20, lmm20 / (disc * C * notional) );
    print("Black implied vol of simulation is ", BlackImpVol/sqrt(Ta) );
    
    return 0;
}
