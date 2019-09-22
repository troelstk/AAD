//
//  main.cpp
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//


#include <iostream>
#include "armadillo.hpp"

//using namespace arma;

#include <vector>
#include <ctime>



#include "AAD.h"
#include "Heston.h"
#include "lsm.h"
#include "vasicek.h"
#include "swap.h"

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

    
    arma::arma_rng::set_seed_random();
    // Create a 4x4 random matrix and print it on the screen
    arma::Mat<double> A = arma::randu(4,4);
    std::cout << "A:\n" << A << "\n";
    // Multiply A with his transpose:
    std::cout << "A * A.t() =\n";
    std::cout <<  A * A.t() << "\n";

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
