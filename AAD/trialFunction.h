//
//  trialFunction.h
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef trialFunction_h
#define trialFunction_h

using namespace std;

// Define function
template<class T> T trialFunction (T x[3]) {
    
    auto temp = 2.0*x[0] + x[1]*x[2] + x[2];
    
    for( int i = 1; i<30; ++i){ // max 30000
        temp = temp + x[0] + x[2] * x[1] + 5.12312*x[0]*x[1]*x[2];
    }
    
    return temp;
};


#endif /* trialFunction_h */
