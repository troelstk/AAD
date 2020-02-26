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


using namespace arma;

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
