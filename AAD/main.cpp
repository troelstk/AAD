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

#include "number.h"
#include "trialFunction.h"
#include "blockList.h"
#include "memoryConsumption.h"

using namespace std;

int main(int argc, const char * argv[]) {
    size_t currentSize = getCurrentRSS( );
    cout << "Current size is " << currentSize << " kb" << endl;
    number x[3] = {1.1, 1.2, 1.3};
    
    number y = trialFunction(x);
    

 
    
 
    y.propagateAdjoints();
    
    for( int i = 0; i<3; i++) {
        cout << "Differential is " << x[i].adjoint() <<  endl;
    }
    

    
    
    // 2.0*x[0] + x[2]*x[1] + x[2] + x[0]
    
    currentSize = getCurrentRSS( );
    cout << "Current size is " << currentSize << " bytes" << endl;
    
    size_t peakSize    = getPeakRSS();
    cout << peakSize << endl;
    
    // << "\n";
    return 0;
}
