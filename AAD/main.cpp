//
//  main.cpp
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//


#include <iostream>
//#include "armadillo.hpp"

//using namespace arma;

#include <vector>
#include <ctime>



#include "AAD.h"
//#include "trialFunction.h"
#include "Heston.h"
#include "lsm.h"
#include "vasicek.h"
#include "swap.h"

using namespace std;



int main(int argc, const char * argv[]) {
   
    //test_heston();
    
    /*clock_t begin_time = clock();
    lsm(200000, 50);
    auto time_aad =  float( clock () - begin_time )/  CLOCKS_PER_SEC;
    cout << "Calculation time: " << time_aad  << " seconds" << endl;
    
    
    double res = bsMonteCarlo(100, 0.2, 1, 100, 123, 234);
    cout << "Result is " <<  res <<  endl;*/

    double mat = 4.0;
    int seed1 = 12, seed2 = 2344, nPaths = 10000;
    int nSteps = 25;
    double r0 = 0.025, t = 0.0, yearly = 1.0 ;
    
    double k = 0.2; // Speed of mean reversion
    double theta = 0.025; // Long term average rate
    double sig = 0.02; // instantaneous volatility
    
    vector<double> params = {k, theta, sig};
    
    double swap_rate_1 = swap_rate(r0, t, t, mat, params, P, yearly);
    cout << "Swap rate is " << swap_rate_1 << endl;
    
    double r_fix = 0.025, r_short = 0.02, notional = 100, Ta = 1, Tb = 4, yearly_payments = 1;
    
    double swap_price1 = swap_price(r_fix, r_short, notional, t, Ta, Tb, params, P, yearly_payments);
    
    cout << "Swap price is " << swap_price1 << endl;
    
    //vasicek_swap(vector<T> params, T t, T Ta, T Tb, T r_fix, T notional, T r0,
    //             int seed1, int seed2, int nPaths, int nSteps, T yearly_payments)
    double sim_swap = vasicek_swap(params, t, Ta, Tb, r_fix, notional, r0,
                                   seed1, seed2, nPaths, nSteps, yearly_payments);
    cout << "simulated swap value is " << sim_swap << endl;
    
    
    
    
    /*cout << "Yield curve: " << endl;
    
    params[0] = 4;
    for(int i=1; i<4*50; ++i){
        double DF = P(0.01, 0.0, double(i)*0.25, params);
        double yield = - log(DF)/(double(i)*0.25-0.0);
        cout << yield << endl;
    }*/
    
    
    /*swap(T r, T t, T Ta, T Tb,
         vector<T> params,
         T (*P)(T r_, T t_, T T_, vector<T> params_),  int nSteps)*/
    
    //double res = vasicek(k, theta, mat, sig, seed1, seed2, nPaths, nSteps);
    //cout << "Result is " <<  res <<  endl;
    
    // << "\n";
    return 0;
}

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
