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
    int nSteps = 1, nPaths = 5;
    double t = 0.0;
    
    double Ta = 1.0, Tb = 4.0, yearly_payments = 1.0;
    double strike = 0.0, notional = 100, r_fix = 0.025;
    
    int M = 3; // Dimension of forward curve
    vector<vector<double>> vol(M,  vector<double>(M));
    vector<vector<double>> corr(M, vector<double>(M));
    vector<double> F = {0.02, 0.025, 0.03};
    
    print( "Swap rate is ", SR_from_F(F, yearly_payments, Ta, Tb) );
    
    // Table 3: Constant vol for each F regardless of time t
    vol[0] = {0.20, 0.00, 0.00};
    vol[1] = {0.25, 0.20, 0.00};
    vol[2] = {0.30, 0.30, 0.20};
    // Burde bruge table 5 med time dependent vol
    
    corr[0] = {1.0, 0.8, 0.7};
    corr[1] = {0.8, 1.0, 0.8};
    corr[2] = {0.7, 0.8, 1.0};
    // Bør lave mere sofistikeret corr
    
    double lmm = LMM_swaption(vol, corr, F, t, Ta, Tb, r_fix, notional,
                          seed1, seed2, nPaths, nSteps, yearly_payments, strike);
    cout << "Swaption price is " << lmm << endl;
    
    /*vec means(5, fill::zeros);
    mat B(5, 5, fill::eye);
    mat C = B.t() * B;
    
    print("Cov is ", C);
    mat X = arma::mvnrnd(means, C, 2);
     */

    return 0;
}

/*cout << "Yield curve: " << endl;
 params[0] = 2;
 params[1] = 0.04;
 for(int i=1; i<4*20; ++i){
 double DF = P(0.01, 0.0, double(i)*0.25, params);
 double yield = yield_from_df(DF, 1.0, 0.0, double(i) * 0.25);
 cout << yield << endl;
 }*/



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
