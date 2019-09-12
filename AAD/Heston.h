#pragma once

#include <math.h>
#include "gaussians.h"
#include <cstdlib>
#include "AAD.h"

template<class T>
T Heston_MC_AAD(T k, T lambda, T eps, T rho, T spot, T K, T TimeToMat,
                double nPaths, T nSteps, int seed = 1) {
    //int int_nSteps = (int)(nSteps);
    int int_nPaths = (int)(nPaths);
    T dt = TimeToMat / nSteps;
    if (seed != 1) {
        srand((int)(seed));
    } 
    
    double N;
    T v, S, res(0.0), payoff(0.0), a = rho, b = sqrt(1 - rho * rho), x, y;
    
    number::tape->mark();
    
    for (int n = 0; n < int_nPaths; n++) {
        number::tape->rewindToMark();
        S = spot;
        v = 1;
        
        for (int i = 0; i < nSteps; i++) { //
            double unif = double(rand() + 1.0) / double(RAND_MAX + 2.0);
            N = invNormalCdf(unif);
            number N2 = a * N + b * invNormalCdf(double(rand() + 1.0) / double(RAND_MAX + 2.0));
            x = 1 + exp(-k * dt) * (v - 1);
            y = sqrt(log(1 + eps * eps * ( v * (exp(- k * dt) - exp(- 2 * k * dt))
                                          + 0.5 * pow((1 - exp(-k * dt)), 2.0))
                                          / (k * x * x)));
            S += lambda * sqrt(v) * S * sqrt(dt) * N2;
            v = x * exp(- y * y * 0.5 + y * N);
        }
        number payoff(0.0);
        if (S > K) {
            payoff = S - K;
        }
        payoff.propagateToMark();
        res += payoff;
    }
    //number::propagateMarkToStart();
    
    //res.propagateMarkToStart();
    
    return res / nPaths;
}

template<class T>
T Heston_MC(T k, T lambda, T eps, T rho, T spot, T K, T TimeToMat, double nPaths, T nSteps, int seed = 1) {
    //int int_nSteps = (int)(nSteps);
    int int_nPaths = (int)(nPaths);
    T dt = TimeToMat / nSteps;
    if (seed != 1) {
        srand((int)(seed));
    }
    
    double N;
    T v, S, res(0.0), payoff(0.0), a = rho, b = sqrt(1 - rho * rho), x, y;
    for (int n = 0; n < int_nPaths; n++) {
        
        S = spot;
        v = 1;
        for (int i = 0; i < nSteps; i++) {
            double unif = double(rand() + 1.0) / double(RAND_MAX + 2.0);
            N = invNormalCdf(unif);
            double N2 = a * N + b * invNormalCdf(double(rand() + 1.0) / double(RAND_MAX + 2.0));
            x = 1 + exp(-k * dt)*(v - 1);
            y = sqrt(log(1 +
                         eps * eps*(v*(exp(-k * dt) - exp(-2 * k * dt)) + 0.5*pow((1 - exp(-k * dt)), 2.0))
                         / (k*x*x)
                         ));
            S += lambda * sqrt(v)*S*sqrt(dt)*N2;
            v = x * exp(-y * y*0.5 + y * N);
        }
        
        if (S > K) {
            payoff = S - K;
        }
        else {
            payoff = 0.0;
        }
        res += payoff;
    }
    
    return res / nPaths;
}
// Define function
template<class T> T trialFunction (T x[3]) {
    
    auto temp = x[0] + x[1]*x[2] + x[2];
    
    for( int i = 1; i<1000; ++i){ // max 30000
        temp += x[0] + x[2] * x[1] + 5.12312 + x[0] *x[1] *x[2];
    }
    
    return temp;
};

void test_heston() {
    
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
     cout << "adj is " << k.adjoint() << endl;
     */
    
    clock_t begin_time = clock();
    
    number::tape->rewind();
    
    double k1 = 1, lambda1 = 0.2, eps1 = 0.5, rho1 = 0.0, spot1=100, strike1 = 100, timeToMat1 = 1.0;
    double nPaths1 = 10000, nSteps1 = 20;
    
    number heston_res;
    number k(k1), lambda(lambda1), eps(eps1), rho(rho1), spot(spot1), strike(strike1), timeToMat(timeToMat1), nSteps(nSteps1);
    double nPaths = nPaths1;
    
    int seed = 123;
    
    begin_time = clock();
    heston_res = Heston_MC_AAD(k, lambda, eps, rho, spot, strike, timeToMat, nPaths, nSteps, seed);
    cout << "Heston result with AAD: " << heston_res.value() << endl;
    auto time_aad =  float( clock () - begin_time )/  CLOCKS_PER_SEC;
    cout << "Calculation time: " << time_aad  << " seconds" << endl;
    
    
    cout << "Delta is " << spot.adjoint()/double(nPaths) << endl;
    
    // Heston without AAD:
    double heston_res1;
    begin_time = clock();
    heston_res1 = Heston_MC(k1, lambda1, eps1, rho1, spot1, strike1, timeToMat1, nPaths1, nSteps1, seed);
    cout << "Heston result without AAD: " << heston_res1 << endl;
    auto time_without_aad = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    cout << "Calculation time: " << time_without_aad << " seconds" << endl;
    
    cout << "AAD is " << time_aad/time_without_aad << " times slower\n";
    
    
    double bump = 0.0001;
    // Central difference delta:
    double heston_res2 =
    Heston_MC(k1, lambda1, eps1, rho1, spot1+bump, strike1, timeToMat1, nPaths1, nSteps1, seed);
    cout << "Delta should be (FD): " << (heston_res2-heston_res1)/(bump) << endl;
    
    
}
