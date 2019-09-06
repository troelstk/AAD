//
//  main.cpp
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//


#include <iostream>

#include <vector>
#include <ctime>



#include "AAD.h"
//#include "trialFunction.h"
#include "Heston.h"

using namespace std;

// Define function
template<class T> T trialFunction (T x[3]) {
    
    auto temp = x[0] + x[1]*x[2] + x[2];
    
    for( int i = 1; i<1000; ++i){ // max 30000
        temp += x[0] + x[2] * x[1] + 5.12312 + x[0] *x[1] *x[2];
    }
    
    return temp;
};

int main(int argc, const char * argv[]) {
    size_t currentSize = getCurrentRSS( );
    cout << "Current size is " << currentSize/10e6 << " MB" << endl;
    
    number::tape->rewind();
    
    number x[3] = {number(1.1), number(1.2), number(1.3)};
    clock_t begin_time = clock();
    
    number y = trialFunction(x);
    cout << "Result for number's: " << y.value() << endl;
    y.propagateToStart();
    cout << "Propagation time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    for( int i = 0; i<3; i++) {
        cout << "Differential is " << x[i].adjoint() <<  endl;
    }
    size_t diff =  getCurrentRSS( )-currentSize;
    cout << "Additional memory used " << diff/10e6 << " MB" << endl;
    
    
    double x_double[3] = {1.1, 1.2, 1.3};
    begin_time = clock();
    double y_double = trialFunction(x_double);
    cout << "Result for doubles: " << y_double << endl;
    cout << "Calculation time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    
    
    
    
    double heston_res;
    double k = 0.15, lambda = 0.2, eps = 0.5, rho = 0.5, spot=100, strike = 100, timeToMat = 1;
    double nPaths = 100, nSteps = 20;
    heston_res = Heston_MC(k, lambda, eps, rho, spot, strike, timeToMat, nPaths, nSteps, 123);
    
    cout << heston_res << endl;
    
    //number x = number(3.09);
    
    number::tape->rewind();
    
    // << "\n";
    return 0;
}
