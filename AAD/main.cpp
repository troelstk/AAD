//
//  main.cpp
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//


#include "/usr/local/Cellar/armadillo/9.700.2/include/armadillo"

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

    
    int seed1 = 132, seed2 = 1222;
    int nSteps = 50, nPaths = 10;
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

    // LMM_swaption(vector<vector<T>> & vol, vector<vector<T>> & corr, vector<T> & initF, T t, T Ta, T Tb, T r_fix, T notional, int seed1, int seed2, int nPaths, int nSteps, T yearly_payments)
    double lmm = LMM_swaption(vol, corr, F, t, Ta, Tb, r_fix, notional,
                          seed1, seed2, nPaths, nSteps, yearly_payments);
    auto time =  float( clock () - begin_time )/ CLOCKS_PER_SEC;
    print("LMM Swaption price is ", lmm);
    print("Calculation time: ", time, " seconds");

    
    double sum = 0.0;
    for(int i = 0; i<F.size(); ++i){
        double weight = w(F, 1.0, double(i), int(Ta));
        //print("weight is ", weight);
        sum += weight * F[i];
    }
    print("swap rate is ", sum);
    
    double rebonato = sqrt(vol_TFM(F, yearly_payments, Ta, corr, vol, swap_rate, Ta) );
    print("Rebonato vol is ", rebonato);
    
    double Cab = C_ab( F, yearly_payments);
    
    // Black(T K, T F0, T vol)
    double black = Cab * BlackCall( r_fix, swap_rate, rebonato);

    print("Black price ", black*100);

    
    
    
    
    
    
    // Test case as in 8.2 in Brigo
    
    print("TEST CASE: ");
    vector<double> F20 = {0.0469, 0.0501, 0.0560, 0.0584, 0.0600, 0.0613, 0.0628, 0.0627, 0.0629, 0.0623,
                          0.0630, 0.0636, 0.0643, 0.0648, 0.0653, 0.0640, 0.0630, 0.0618, 0.0607, 0.0594 };
    
    // Test 1.b, formulation 7 with a=0, b=0, c=1, d=0, all correlations set to 1
    M = 20;
    vector<vector<double>> vol20(M,  vector<double>(M));
    vector<vector<double>> corr20(M, vector<double>(M));

    double swap_rate20 = SR_from_F(F20, yearly_payments, 10, 20);
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

    //print(vol20);
    //print(corr20);
    int start_idx = 10;
    Ta = 10;
    Tb = 20;
    double rebonato2 = sqrt(vol_TFM(F20, yearly_payments, Ta, corr20, vol20, swap_rate20, start_idx) );
    print("Rebonato vol is ", rebonato2/sqrt(Ta));
    
    double C = C_ab( F20, yearly_payments);
    
    // Black(T K, T F0, T vol)
    double black20 = C * BlackCall( swap_rate20, swap_rate20, rebonato2);
    
    print("Black price ", black20 * 100);
    
    nSteps = 50;
    //nPaths = 10000;
    
    double lmm20 = LMM_swaption(vol20, corr20, F20, t, Ta, Tb, swap_rate20, notional,
                          seed1, seed2, nPaths, nSteps, yearly_payments);
    time =  float( clock () - begin_time )/ CLOCKS_PER_SEC;
    print("LMM Swaption price is ", lmm20);
    print("Calculation time: ", time, " seconds");
                            
    
    return 0;
}

/*cout << "Yield curve: " << endl;
 params[0] = 2;
 params[1] = 0.04;
 for(int i=1; i<4*20; ++i){
 double DF = P(0.01, 0.0, double(i)*0.25, params);
 double yield = yield_from_df(DF, 1.0, 0.0, double(i) * 0.25);
 cout << yield << endl;
 }
 double weight = w(F, 1.0, 1.0);
 print("Weight is ", weight);
 weight += w(F, 1.0, 2.0);
 print("Weight 2 is ", weight);
 weight += w(F, 1.0, 3.0);
 print("Weight 3 is ", weight);
 */



/*arma::arma_rng::set_seed_random();
// Create a 4x4 random matrix and print it on the screen
arma::Mat<double> A = arma::randu(4,4);
std::cout << "A:\n" << A << "\n";
// Multiply A with his transpose:
std::cout << "A * A.t() =\n";
std::cout <<  A * A.t() << "\n";*/

/*#include <iostream>

int main(int argc, const char **argv) {
    // Initialize the random generator
    arma::arma_rng::set_seed_random();
    // Create a 4x4 random matrix and print it on the screen
    arma::Mat<double> A = arma::randu(4,4);
    std::cout << "A:\n" << A << "\n";
    // Multiply A with his transpose:
    std::cout << "A * A.t() =\n";
    std::cout << A * A.t() << "\n";
    // Access/Modify rows and columns from the array:
    A.row(0) = A.row(1) + A.row(3);
    A.col(3).zeros();
    std::cout << "add rows 1 and 3, store result in row 0, also fill 4th column with zeros:\n";
    std::cout << "A:\n" << A << "\n";
    // Create a new diagonal matrix using the main diagonal of A:
    arma::Mat<double>B = arma::diagmat(A);
    std::cout << "B:\n" << B << "\n";
    // Save matrices A and B:
    A.save("A_mat.txt", arma::arma_ascii);
    B.save("B_mat.txt", arma::arma_ascii);
    return 0;
}*/
