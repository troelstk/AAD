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
    //size_t currentSize = getCurrentRSS( );
    //cout << "Current size is " << currentSize/10e6 << " MB" << endl;
    /*number x[3] = {number(1.1), number(1.2), number(1.3)};
    
    
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
    cout << "Calculation time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;*/
    
    /*number result(0.0), temp(0.0);
     
     number::tape->mark();
     for(int i=0; i<5; ++i){
     number::tape->rewindToMark();
     number temp = lambda * k;
     temp.propagateToMark();
     result += temp;
     cout << i << endl;
     }
     
     cout << "val is " << result.value() << endl;
     cout << "adj is " << k.adjoint() << endl;*/
    
    clock_t begin_time = clock();
    
    number::tape->rewind();
    
    double k1 = 1, lambda1 = 0.2, eps1 = 0.5, rho1 = 0.0, spot1=100, strike1 = 100, timeToMat1 = 1.0;
    double nPaths1 = 100000, nSteps1 = 20;
    
    number heston_res;
    number k(k1), lambda(lambda1), eps(eps1), rho(rho1), spot(spot1), strike(strike1), timeToMat(timeToMat1), nSteps(nSteps1);
    double nPaths = nPaths1;

    
    begin_time = clock();
    heston_res = Heston_MC_AAD(k, lambda, eps, rho, spot, strike, timeToMat, nPaths, nSteps, 12);
    cout << "Heston result with AAD: " << heston_res.value() << endl;
    cout << "Calculation time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds" << endl;
    
    
    cout << "Delta is " << spot.adjoint()/double(nPaths) << endl;
    
    // Heston without AAD:
    double heston_res1;
    begin_time = clock();
    heston_res1 = Heston_MC(k1, lambda1, eps1, rho1, spot1, strike1, timeToMat1, nPaths1, nSteps1, 12);
    cout << "Heston result without AAD: " << heston_res1 << endl;
    cout << "Calculation time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds" << endl;
    
    
    
    // Central difference delta:
    double heston_res2 =
        Heston_MC(k1, lambda1, eps1, rho1, spot1+0.001, strike1, timeToMat1, nPaths1, nSteps1, 12);
    cout << "Delta should be (CD): " << (heston_res2-heston_res1)/(0.001) << endl;
    
    
    
    
    // << "\n";
    return 0;
}
