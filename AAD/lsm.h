//
//  lsm.h
//  AAD
//
//  Created by Troels Tang Karlsen on 12/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef lsm_h
#define lsm_h
#include <iostream>
#include <random>
#include <cstdlib>

#include "gaussians.h"

using namespace std;

#include "mrg32k3a.h"
// #include "armadillo.hpp"


using namespace arma;

double lsm(double nPaths, double nSteps ) {
    
    double T = 1.0, vol = 0.2, dt=T/nSteps, spot = 100, strike = 100;
    
    

    //double dt = mat / double(nSteps);
    double sdt = vol * sqrt(dt);
    double halfvar = 0.5*sdt*sdt;
    double logS0 = log(spot);
    double res = 0.0;
    
    Mat<double> S = mat(int(nPaths), int(nSteps), fill::none);
    
    S.col(0).fill(logS0);
    
    /*mrg32k3a myRNG ;
    myRNG.init(nSteps);
    
    vector<double> gauss(nSteps);*/
    
    // Fill out matrix
    for (int n = 0; n < nPaths; n++) {
        for (int i = 1; i < nSteps; i++) {
            double unif = double(rand() + 1.0) / double(RAND_MAX + 2.0);
            double N = invNormalCdf(unif);
            //logS += - halfvar + sdt * N;
            S(n,i) = S(n,i-1) - halfvar + sdt * N;
        }
        //cout << exp(S(n,nSteps-1)) << endl;
    }
    Col<double> strikes(nPaths);
    strikes.fill(strike);
    
    Col<double> scol  = S.col(nSteps-1);
    umat res2 = ( scol >= strikes );
    
    

    /*spotatT = exp(logS);
    if (spotatT > strike) {
        payoff = spotatT - strike;
    }
    else {
        payoff = 0.0;
    }
    res += payoff;*/
    
    return res/double(nPaths);
}


double bsMonteCarlo(double spot, double vol, double mat, double strike, int seed1, int seed2)
{
    //srand(23453);
    
    int nSteps = 20;
    int nPaths = 10000;
    double dt = mat / double(nSteps);
    double sdt = vol * sqrt(dt);
    double halfvar = 0.5 * sdt * sdt;
    double logS0 = log(spot);
    double res = 0.0;
    double spotatT = 0.0;
    double payoff = 0.0;
    
    mrg32k3a myRNG(seed1, seed2) ;
    myRNG.init(nSteps);

    vector<double> gauss(nSteps);
    
    for (int n = 0; n < nPaths; n++) {
        double logS = logS0;
        myRNG.nextG(gauss);
        
        for (int i = 0; i < nSteps; i++) {
            //double unif = double(rand() + 1.0) / double(RAND_MAX + 2.0);
            double N = invNormalCdf(gauss[i]);
            logS += - halfvar + sdt * N;
        }
        spotatT = exp(logS);
        if (spotatT > strike) {
            payoff = spotatT - strike;
        }
        else {
            payoff = 0.0;
        }
        res += payoff;
    }
    
    return res/double(nPaths);
}


#endif /* lsm_h */
