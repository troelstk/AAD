//
//  node.h
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef node_h
#define node_h

using namespace std;

#include <string>
#include <vector>
#include <memory>

class node {
    
    friend class Tape;
    friend class number;
    
    // Number of childs
    const size_t n;
    
    // Adjoint of this node
    double mAdjoint = 0;
    
    // Pointers to derivatives between child nodes and this node
    double* pDerivatives;
    
    // Pointers to adjoints of child nodes 
    double** pAdjPtrs;

public:
    // Initialize with n
    node(const size_t N = 0) : n(N) { }
    
    double& adjoint() {
        return mAdjoint;
    }
    
    // Method to add to node w_j's childrens adjoints
    void propagateOne()
    {
        // Return if zero childs (nothing to propagate) or mAdjoint is 0 (no need to propagate as all terms added are just 0)
        if ( !n || !mAdjoint) return;
        
        for(size_t i = 0; i < n; ++i ) {
            *(pAdjPtrs[i]) += pDerivatives[i] * mAdjoint;
        }
    }
};



#endif /* node_h */
