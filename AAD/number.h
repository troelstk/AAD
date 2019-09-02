//
//  number.h
//  DAG
//
//  Created by Troels Tang Karlsen on 29/08/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef number_h
#define number_h

#include "node.h"

#include<memory>



class number {

    // a number holds a pointer to a node
    node* myNode;
    
public:
    // The tape: A vector of pointers to nodes
    static vector<unique_ptr<node>> tape;
    
    // Create new leaf to hold the number
    number(double val) : myNode(new leaf(val))
    {
        // Put this numbers node on the tape
        tape.push_back(unique_ptr<node>(myNode));
    };
    
    // Construct a number from a shared pointer to a node,
    number(node * inp_node) : myNode(inp_node) {}
    
    // Method to return a numbers node
    node* getNode() {
        return myNode;
    }
    
    // Set value of a number, only leaves can hold a value
    void setVal(double val){
        dynamic_cast<leaf*>(myNode)->setValue(val) ;
        // returns a copy of myNode with its stored pointer casted dynamically from node to leaf
    }
    
    // Get value, only from leaves
    double getVal(){
        return dynamic_cast<leaf*>(myNode)->getValue() ;
    }
    
    double& adjoint() {
        //cout << "adjoint is " << myNode->getAdjoint() <<"\n";
        return myNode->adjoint();
    }
    
    void propagateAdjoints()
    {
        // Set all adjoints to 0
        myNode->resetAdjoints();
        
        // Set adjoint of this last node to 1
        myNode->adjoint() = 1.0;
        
        // Start from end of tape and iterate backwards until this node is reached (rbegin points to last)
        auto iter = tape.rbegin();
        while(iter->get() != myNode){
            ++iter;
        }
        // Now iter points to this node on the tape
        while(iter != tape.rend()) {
            (*iter)->propagateAdjoint();
            ++iter; // Jump one backwards
        }
    }
};

vector<unique_ptr<node>> number::tape;

number operator+(number left, number right){
    // Create plus-node
    node* n = new plusNode(left.getNode(), right.getNode());
    
    // Push the new node onto the tape
    number::tape.push_back(unique_ptr<node>(n));
    
    return n;
}

number operator*(number left, number right){
    // Create times-node
    node* n = new timesNode(left.getNode(), right.getNode());
    
    // Push the new node onto the tape
    number::tape.push_back(unique_ptr<node>(n));
    return n;
}

#endif /* number_h */
