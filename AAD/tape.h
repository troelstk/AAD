//
//  tape.h
//  
//
//  Created by Troels Tang Karlsen on 02/09/2019.
//

#ifndef tape_h
#define tape_h

#include "blockList.h"
#include "node.h"
#include <iostream>

// Number of nodes in a block
constexpr size_t blockSize = 16384;
constexpr size_t dataSize = 65536;


class Tape {
    // Blocklist for nodes
    blockList<node, blockSize> myNodes;
    
    // blocklist for derivatives and child adjoint pointers
    blockList<double, dataSize> myDerivs;
    blockList<double*, dataSize> myArgPtrs;
    
    //char myPad[64];
    
public:
    // Put node on blocklist for nodes, derivatives on blocklist for derivs, adjoints on blocklist for adj.. 
    template <size_t N>
    node* recordNode()
    {
        node* Node = myNodes.emplace_back(N);
        //std::cout << "Debug 1 " << endl;
        if constexpr(N != 0)
        {
            // Set derivatives
            Node->pDerivatives = myDerivs.emplace_back_multi<N>();
            // Set child adjoints
            Node->pAdjPtrs = myArgPtrs.emplace_back_multi<N>();
        }
        return Node;
    }
    
    void resetAdjoints()
    {
        for (node& Node : myNodes){
            Node.mAdjoint = 0.0;
        }
    }
    
    void clear()
    {
        myDerivs.clear();
        myArgPtrs.clear();
        myNodes.clear();
    }
    
    void rewind() {
        myDerivs.rewind();
        myArgPtrs.rewind();
        myNodes.rewind();
    }
    
    void mark(){
        myDerivs.setMark();
        myArgPtrs.setMark();
        myNodes.setMark();
    }
    
    void rewindToMark(){
        myDerivs.rewindToMark();
        myArgPtrs.rewindToMark();
        myNodes.rewindToMark();
    }
    
    // iterators
    
    using iterator = blockList<node, blockSize>::iterator;
    
    auto begin()
    {
        return myNodes.begin();
    }
    auto end()
    {
        return myNodes.end();
    }
    auto markIt()
    {
        return myNodes.mark();
    }
    auto find(node* Node)
    {
        return myNodes.find(Node);
    }
    
};


#endif /* tape_h */
