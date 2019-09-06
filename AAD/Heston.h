#pragma once

#include <math.h>
#include "gaussians.h"

template<class T>
T Heston_MC(T k, T lambda, T eps, T rho, T spot, T K, T TimeToMat, double nPaths, double nSteps, double seed = 1) {
    int int_nSteps = (int)(nSteps);
    int int_nPaths = (int)(nPaths);
    double dt = TimeToMat / nSteps;
    if (seed != 1) {
        srand((int)(seed));
    }
    
    double N, N2;
    T v, S, res = 0.0, payoff = 0.0, a = rho, b = sqrt(1 - rho * rho), x, y;
    for (int n = 0; n < int_nPaths; n++) {
        S = spot;
        v = 1;
        for (int i = 0; i < int_nSteps; i++) {
            N = invNormalCdf(double(rand() + 1) / (RAND_MAX + 2));
            N2 = a * N + b * invNormalCdf(double(rand() + 1) / (RAND_MAX + 2));
            x = 1 + exp(-k * dt)*(v - 1);
            y = sqrt(log(1 +
                         eps * eps*(v*(exp(-k * dt) - exp(-2 * k * dt)) + 0.5*pow((1 - exp(-k * dt)), 2))
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
