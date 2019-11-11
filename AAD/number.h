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
#include "tape.h"
#include "gaussians.h"

#include <memory>
#include <iostream>
#include <cmath>

using namespace std;


class number {

    // a number holds a pointer to a node
    node* myNode;
    
    double myValue;
    
    //// Access to Derivatives and adjoints of this numbers node
    // For unary functions
    double& derivative() {
        return myNode->pDerivatives[0];
    }
    // For binary functions
    double& lDer() {
        return myNode->pDerivatives[0];
    }
    double& rDer() {
        return myNode->pDerivatives[1];
    }
    
    
    // For unary functions
    double*& adj() {
        return myNode->pAdjPtrs[0];
    }
    // For binary functions
    double*& lAdj() {
        return myNode->pAdjPtrs[0];
    }
    double*& rAdj() {
        return myNode->pAdjPtrs[1];
    }
    
public:
    
    static thread_local Tape* tape;
    
    double& value(){
        return myValue;
    }
    double value() const
    {
        return myValue;
    }
    
    // Constructors:
    number() {}
    
    explicit number (const double val): myValue(val) {
        createNode<0>();
    }
    
    // Create this numbers node on the tape
    template<size_t N>
    void createNode()
    {
        myNode = tape->recordNode<N>();
    }
    
    // Assignment of numeric value creates a leaf node
    number& operator=(const double & val){
        myValue = val;
        createNode<0>();
        return *this;
    }
    
    // Put a number on tape explicitly (eg. after a rewind/clear of the tape)
    void putOnTape() {
        createNode<0>();
    }
    
    // Method to return a numbers node
    node* getNode() {
        return myNode;
    }
    
    
    double& adjoint() {
        //cout << "adjoint is " << myNode->getAdjoint() <<"\n";
        return myNode->adjoint();
    }
    
    explicit operator double& () {return myValue;};
    explicit operator double() const {return myValue;};
    
    static void propagateAdjoints(Tape::iterator propagateFrom,
                                  Tape::iterator propagateTo) {
        // Propagate adjoints from back to start
        auto it = propagateFrom;
        while(it != propagateTo) {
            it->propagateOne();
            --it;
        }
        it->propagateOne();
    }
    
    // Overload: Set adjoint to 1 and propagate back to propagateTo
    void propagateAdjoints(Tape::iterator propagateTo) {
        // Propagate adjoints from back to start
        adjoint() = 1.0;
        auto propagateFrom = tape->find(myNode);
        
        propagateAdjoints(propagateFrom, propagateTo);
    }
    void propagateToStart() {
        // Propagate adjoints to start
        propagateAdjoints(tape->begin());
    }
    void propagateToMark() {
        // Propagate adjoints to mark
        propagateAdjoints(tape->markIt());
    }

    static void propagateMarkToStart() {
        // Propagate adjoints from mark-1 to start
        propagateAdjoints(prev(tape->markIt()), tape->begin());
    }
    
    // Operator overloads
    
    
    /* The following operator overloads are strict IP of Antoine Savine */
    
    inline friend number operator*(const number& lhs, const number& rhs)
    {
        const double e = lhs.value() * rhs.value();
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), rhs.Node(), e);
        //  Eagerly compute derivatives
        result.lDer() = rhs.value();
        result.rDer() = lhs.value();
        
        return result;
    }
    inline friend number operator*(const number& lhs, const double& rhs)
    {
        const double e = lhs.value() * rhs;
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = rhs;
        
        return result;
        
    }
    inline friend number operator*(const double& lhs, const number& rhs)
    {
        return rhs * lhs;
    }
    
    inline friend number operator+(const number& lhs, const number& rhs)
    {
        const double e = lhs.value() + rhs.value();
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), rhs.Node(), e);
        //  Eagerly compute derivatives
        result.lDer() = 1.0;
        result.rDer() = 1.0;
        
        return result;
    }
    inline friend number operator+(const number& lhs, const double& rhs)
    {
        const double e = lhs.value() + rhs;
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = 1.0;
        
        return result;
        
    }
    inline friend number operator+(const double& lhs, const number& rhs)
    {
        return rhs + lhs;
    }
    
    inline friend number operator-(const number& lhs, const number& rhs)
    {
        const double e = lhs.value() - rhs.value();
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), rhs.Node(), e);
        //  Eagerly compute derivatives
        result.lDer() = 1.0;
        result.rDer() = -1.0;
        
        return result;
    }
    inline friend number operator-(const number& lhs, const double& rhs)
    {
        const double e = lhs.value() - rhs;
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = 1.0;
        
        return result;
        
    }
    inline friend number operator-(const double& lhs, const number& rhs)
    {
        const double e = lhs - rhs.value();
        //  Eagerly evaluate and put on tape
        number result(rhs.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = -1.0;
        
        return result;
    }
    inline friend number operator/(const number& lhs, const number& rhs)
    {
        const double e = lhs.value() / rhs.value();
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), rhs.Node(), e);
        //  Eagerly compute derivatives
        const double invRhs = 1.0 / rhs.value();
        result.lDer() = invRhs;
        result.rDer() = -lhs.value() * invRhs * invRhs;
        
        return result;
    }
    inline friend number operator/(const number& lhs, const double& rhs)
    {
        const double e = lhs.value() / rhs;
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = 1.0 / rhs;
        
        return result;
        
    }
    inline friend number operator/(const double& lhs, const number& rhs)
    {
        const double e = lhs / rhs.value();
        //  Eagerly evaluate and put on tape
        number result(rhs.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = -lhs / rhs.value() / rhs.value();
        
        return result;
    }
    
    inline friend number pow(const number& lhs, const number& rhs)
    {
        const double e = pow(lhs.value(), rhs.value());
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), rhs.Node(), e);
        //  Eagerly compute derivatives
        result.lDer() = rhs.value() * e / lhs.value();
        result.rDer() = log(lhs.value()) * e;
        
        return result;
    }
    inline friend number pow(const number& lhs, const double& rhs)
    {
        const double e = pow(lhs.value(), rhs);
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = rhs * e / lhs.value();
        
        return result;
    }
    inline friend number pow(const double& lhs, const number& rhs)
    {
        const double e = pow(lhs, rhs.value());
        //  Eagerly evaluate and put on tape
        number result(rhs.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = log(lhs) * e;
        
        return result;
    }
    
    inline friend number max(const number& lhs, const number& rhs)
    {
        const bool lmax = lhs.value() > rhs.value();
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), rhs.Node(), lmax? lhs.value() : rhs.value());
        //  Eagerly compute derivatives
        if (lmax)
        {
            result.lDer() = 1.0;
            result.rDer() = 0.0;
        }
        else
        {
            result.lDer() = 0.0;
            result.rDer() = 1.0;
        }
        
        return result;
    }
    inline friend number max(const number& lhs, const double& rhs)
    {
        const bool lmax = lhs.value() > rhs;
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), lmax ? lhs.value() : rhs);
        //  Eagerly compute derivatives
        result.derivative() = lmax ? 1.0 : 0.0;
        
        return result;
    }
    inline friend number max(const double& lhs, const number& rhs)
    {
        const bool rmax = rhs.value() > lhs;
        //  Eagerly evaluate and put on tape
        number result(rhs.Node(), rmax ? rhs.value() : lhs);
        //  Eagerly compute derivatives
        result.derivative() = rmax ? 1.0 : 0.0;
        
        return result;
    }
    
    inline friend number min(const number& lhs, const number& rhs)
    {
        const bool lmin = lhs.value() < rhs.value();
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), rhs.Node(), lmin ? lhs.value() : rhs.value());
        //  Eagerly compute derivatives
        if (lmin)
        {
            result.lDer() = 1.0;
            result.rDer() = 0.0;
        }
        else
        {
            result.lDer() = 0.0;
            result.rDer() = 1.0;
        }
        
        return result;
    }
    inline friend number min(const number& lhs, const double& rhs)
    {
        const bool lmin = lhs.value() < rhs;
        //  Eagerly evaluate and put on tape
        number result(lhs.Node(), lmin ? lhs.value() : rhs);
        //  Eagerly compute derivatives
        result.derivative() = lmin ? 1.0 : 0.0;
        
        return result;
    }
    inline friend number min(const double& lhs, const number& rhs)
    {
        const bool rmin = rhs.value() < lhs;
        //  Eagerly evaluate and put on tape
        number result(rhs.Node(), rmin ? rhs.value() : lhs);
        //  Eagerly compute derivatives
        result.derivative() = rmin ? 1.0 : 0.0;
        
        return result;
    }
    
    number& operator+=(const number& arg)
    {
        *this = *this + arg;
        return *this;
    }
    number& operator+=(const double& arg)
    {
        *this = *this + arg;
        return *this;
    }
    
    number& operator-=(const number& arg)
    {
        *this = *this - arg;
        return *this;
    }
    number& operator-=(const double& arg)
    {
        *this = *this - arg;
        return *this;
    }
    
    number& operator*=(const number& arg)
    {
        *this = *this * arg;
        return *this;
    }
    number& operator*=(const double& arg)
    {
        *this = *this * arg;
        return *this;
    }
    
    number& operator/=(const number& arg)
    {
        *this = *this / arg;
        return *this;
    }
    number& operator/=(const double& arg)
    {
        *this = *this / arg;
        return *this;
    }
    
    //  Unary +/-
    number operator-() const
    {
        return 0.0 - *this;
    }
    number operator+() const
    {
        return *this;
    }
    
    //  Overloading continued, unary functions
    
    inline friend number exp(const number& arg)
    {
        const double e = exp(arg.value());
        //  Eagerly evaluate and put on tape
        number result(arg.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = e;
        
        return result;
    }
    
    inline friend number log(const number& arg)
    {
        const double e = log(arg.value());
        //  Eagerly evaluate and put on tape
        number result(arg.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = 1.0 / arg.value();
        
        return result;
    }
    
    inline friend number sqrt(const number& arg)
    {
        const double e = sqrt(arg.value());
        //  Eagerly evaluate and put on tape
        number result(arg.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = 0.5 / e;
        
        return result;
    }
    inline friend number cos(const number& arg)
    {
        const double e = cos(arg.value());
        //  Eagerly evaluate and put on tape
        number result(arg.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = sin(e);
        
        return result;
    }
    
    inline friend number fabs(const number& arg)
    {
        const double e = fabs(arg.value());
        //  Eagerly evaluate and put on tape
        number result(arg.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = arg.value() > 0.0 ? 1.0 : -1.0;
        
        return result;
    }
    
    inline friend number normalDens(const number& arg)
    {
        const double e = normalDens(arg.value());
        //  Eagerly evaluate and put on tape
        number result(arg.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = - arg.value() * e;
        
        return result;
    }
    
    inline friend number normalCdf(const number& arg)
    {
        const double e = normalCdf(arg.value());
        //  Eagerly evaluate and put on tape
        number result(arg.Node(), e);
        //  Eagerly compute derivatives
        result.derivative() = normalDens(arg.value());
        
        return result;
    }
    
    //  Finally, comparison
    
    inline friend bool operator==(const number& lhs, const number& rhs)
    {
        return lhs.value() == rhs.value();
    }
    inline friend bool operator==(const number& lhs, const double& rhs)
    {
        return lhs.value() == rhs;
    }
    inline friend bool operator==(const double& lhs, const number& rhs)
    {
        return lhs == rhs.value();
    }
    
    inline friend bool operator!=(const number& lhs, const number& rhs)
    {
        return lhs.value() != rhs.value();
    }
    inline friend bool operator!=(const number& lhs, const double& rhs)
    {
        return lhs.value() != rhs;
    }
    inline friend bool operator!=(const double& lhs, const number& rhs)
    {
        return lhs != rhs.value();
    }
    
    inline friend bool operator<(const number& lhs, const number& rhs)
    {
        return lhs.value() < rhs.value();
    }
    inline friend bool operator<(const number& lhs, const double& rhs)
    {
        return lhs.value() < rhs;
    }
    inline friend bool operator<(const double& lhs, const number& rhs)
    {
        return lhs < rhs.value();
    }
    
    inline friend bool operator>(const number& lhs, const number& rhs)
    {
        return lhs.value() > rhs.value();
    }
    inline friend bool operator>(const number& lhs, const double& rhs)
    {
        return lhs.value() > rhs;
    }
    inline friend bool operator>(const double& lhs, const number& rhs)
    {
        return lhs > rhs.value();
    }
    
    inline friend bool operator<=(const number& lhs, const number& rhs)
    {
        return lhs.value() <= rhs.value();
    }
    inline friend bool operator<=(const number& lhs, const double& rhs)
    {
        return lhs.value() <= rhs;
    }
    inline friend bool operator<=(const double& lhs, const number& rhs)
    {
        return lhs <= rhs.value();
    }
    
    inline friend bool operator>=(const number& lhs, const number& rhs)
    {
        return lhs.value() >= rhs.value();
    }
    inline friend bool operator>=(const number& lhs, const double& rhs)
    {
        return lhs.value() >= rhs;
    }
    inline friend bool operator>=(const double& lhs, const number& rhs)
    {
        return lhs >= rhs.value();
    }
    
private:
    // Constructors for operator overloads
    number(node& arg, const double val): myValue(val){
        createNode<1>();
        myNode->pAdjPtrs[0] = &arg.mAdjoint;
    }
    
    number(node& lArg, node& rArg, const double val): myValue(val){
        createNode<2>();
        myNode->pAdjPtrs[0] = &lArg.mAdjoint;
        myNode->pAdjPtrs[1] = &rArg.mAdjoint;
    }
    
    // Access this numbers node, by casting away the constness. Needed due to const-incorrectness explained on p. 391
    node& Node() const
    {
        return const_cast<node&>(*myNode);
    }
    
  
};



#endif /* number_h */
