//
//  main.cpp
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//


#include <iostream>
#include "armadillo.hpp"

using namespace arma;

#include <vector>
#include <ctime>



#include "AAD.h"
//#include "trialFunction.h"
#include "Heston.h"
#include "lsm.h"

using namespace std;



int main(int argc, const char * argv[]) {
   
    //test_heston();

    //lsm(50, 20);
    
    double res = bsMonteCarlo(100, 0.2, 1, 100);
    cout << "Result is " <<  res <<  endl;
    
    // << "\n";
    return 0;
}

/*#include <iostream>
#include <armadillo>

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
