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
    
protected:
    // Childs of node
    vector<node*> arguments;
    
    double myResult;
    double myAdjoint = 0.0;

public:
    virtual ~node() {};
    
    double result() {
        return myResult;
    };
    
    double& adjoint() {
        return myAdjoint;
    }
    
    void resetAdjoints() {
        for( auto arg : arguments ) {
            arg->resetAdjoints();
            //cout << "reset adjoint of node with result " << myResult << endl;
        }
        myAdjoint = 0.0;
    }
    
    virtual void propagateAdjoint() = 0;
    
};

class plusNode : public node {
public:
    // Only instructs that a plusNode holds two child nodes
    plusNode(node* left, node* right)
    {
        arguments.resize(2);
        arguments[0] = left;
        arguments[1] = right;
        
        // Eager evalutation
        myResult = left->result() + right->result();
    }
    
    void propagateAdjoint() override
    {
        //cout << "Adjoint is: " << myAdjoint << endl;
        
        arguments[0]->adjoint() += myAdjoint;
        arguments[1]->adjoint() += myAdjoint;
    }
};

class timesNode : public node {
public:
    // Only instructs that a timesNode holds two child nodes
    timesNode(node* left, node* right)
    {
        arguments.resize(2);
        arguments[0] = left;
        arguments[1] = right;
        
        myResult = left->result() * right->result();
    }
    
    void propagateAdjoint() override
    {
        //cout << "Adjoint " << myAdjoint << endl;
        
        arguments[0]->adjoint() += myAdjoint * arguments[1]->result();
        arguments[1]->adjoint() += myAdjoint * arguments[0]->result();
    }
};

class leaf : public node {
    
public:
    
    leaf(double val){
        myResult = val;
    }
    
    void setValue(double val) {
        myResult = val;
    }
    double getValue(){
        return myResult;
    }

    void propagateAdjoint() override { }
    
};



#endif /* node_h */
