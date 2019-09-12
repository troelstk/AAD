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

#include "mrg32k3a.h"
#include "armadillo.hpp"

double lsm(double nPaths, double nSteps ) {
    
    double r = 0.0, q = 0.0, T = 1.0, sigma = 0.2, dt=T/nSteps, BM = 0.0;
    
    Mat<double> S = mat(int(nPaths), int(nSteps), fill::none);
    
    S.col(0).fill(100.0);
    
    //cout << S(0,0) << endl;
    
    for(int i = 0; i<nPaths; ++i){
        //cout << S(i,0) << " " ;
        for(int j = 1; j<nSteps; ++j){
            BM = BM + sqrt(dt)*invNormalCdf(double(rand() + 1.0) / double(RAND_MAX + 2.0));
            S(i,j) = S(i,0) * exp( (r - q - (sigma * sigma) * 0.5) * (j - 1) * dt + sigma * BM);
            //cout << S(i,j) << " ";
        }
        //cout << endl;
    }

    return 0.0;
}


double bsMonteCarlo(double spot, double vol, double mat, double strike)
{
    //srand(23453);
    
    int nSteps = 10;
    int nPaths = 10;
    double dt = mat / double(nSteps);
    double sdt = vol * sqrt(dt);
    double halfvar = 0.5*sdt*sdt;
    double logS0 = log(spot);
    double res = 0.0;
    double spotatT = 0.0;
    double payoff = 0.0;
    
    //mrg32k3a myRNG ;
    //myRNG.init(nSteps);
    
    vector<double> gauss(nSteps);
    
    for (int n = 0; n < nPaths; n++) {
        double logS = logS0;
        //myRNG.nextG(gauss);
        for (int i = 0; i < nSteps; i++) {
            //logS += - halfvar + sdt * invNormalCdf( double(rand() + 1) / (RAND_MAX + 2));
            logS += - halfvar + sdt * gauss[i];
            cout << exp(logS) << endl;
        }
        spotatT = exp(logS);
        cout << spotatT << endl;
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
