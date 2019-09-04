//
//  main.cpp
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//


#include <iostream>

#include <vector>
/*
 #include <ctime>
 const clock_t begin_time = clock();
 // do something
 std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;
*/


#include "AAD.h"
//#include "trialFunction.h"

using namespace std;

// Define function
template<class T> T trialFunction (T x[3]) {
    
    auto temp = x[0] + x[1]*x[2] + x[2];
    
    for( int i = 1; i<2000; ++i){ // max 30000
        temp += x[0] + x[2] * x[1] + 5.12312 + x[0] *x[1] *x[2];
    }
    
    return temp;
};

int main(int argc, const char * argv[]) {
    size_t currentSize = getCurrentRSS( );
    cout << "Current size is " << currentSize << " bytes" << endl;
    
    number::tape->rewind();
    
    number x[3] = {number(1.1), number(1.2), number(1.3)};

    number y = trialFunction(x);
    //y.propagateAdjoints();
    y.propagateToStart();
    for( int i = 0; i<3; i++) {
        cout << "Differential is " << x[i].adjoint() <<  endl;
    }
    size_t diff =  getCurrentRSS( )-currentSize;
    cout << "Additional memory used " << diff << " bytes" << endl;
    
    
    //number x = number(3.09);
    
    number::tape->rewind();
    
    
    // << "\n";
    return 0;
}
