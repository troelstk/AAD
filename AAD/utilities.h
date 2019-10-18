//
//  utilities.h
//  AAD
//
//  Created by Troels Tang Karlsen on 23/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef utilities_h
#define utilities_h
#include <string>

using namespace std;

// Prints the values of a vector in one row
template<class T> void
print(vector<T> input_vec) {
    for( auto &x : input_vec) {
        cout << x << " ";
    }
    cout << endl;
}

// Prints values of a vector of vectors, in one row per vector
template<class T> void
print(vector<vector<T>> input_vec_of_vec) {
    for( auto &x : input_vec_of_vec) {
        print(x);
    }
    cout << endl;
}

// Needed for below print function
struct expand_type {
  template<typename... T>
  expand_type(T&&...) {}
};
// Print variable number of arguments
template<typename... ArgTypes>
void print(ArgTypes... args)
{
    expand_type{ 0, (std::cout << args, 0)... };
    std::cout << endl;
}

// Overload log for vector of T's
template<class T>
vector<T> log(vector<T> input_vec) {
    vector<T> res;
    for( auto &x : input_vec) {
        res.push_back( log(x) );
    }
    return res;
}

// Overload exp for vector of T's
template<class T>
vector<T> exp(vector<T>& input_vec) {
    vector<T> res;
    for( auto &x : input_vec) {
        res.push_back( exp(x) );
    }
    return res;
}
// Sum of a vector
template<class T>
T sum(vector<T>& input_vec) {
    T res(0.0);
    for( auto &x : input_vec) {
        res += x;
    }
    return res;
}

template<class T>
T stdev(vector<T>& input_vec) {
    T stdev(0.0);
    T mean = sum(input_vec)/double(input_vec.size()) ;
    
    for(auto &x : input_vec){
        stdev += pow(x - mean, 2.0) ;
    }
    stdev /= double(input_vec.size());
    
    return sqrt(stdev);
}

template<class T>
T skew(vector<T>& input_vec) {
    T skew(0.0), skewtop(0.0), skewbot(0.0);
    T mean = sum(input_vec)/double(input_vec.size()) ;
    
    for(auto &x : input_vec){
        skewtop +=  pow(x - mean, 3.0)/double(input_vec.size());
        skewbot +=  pow(x - mean, 2.0)/(double(input_vec.size())-1.0);
    }
    skew = skewtop/ pow(skewbot, 1.5);
    
    return skew;
}

template<class T>
T kurtosis(vector<T>& input_vec) {
    T kurtosis(0.0);
    T mean = sum(input_vec)/double(input_vec.size()) ;
    
    for(auto &x : input_vec){
        kurtosis += pow(x - mean, 4.0)/double(input_vec.size());
    }
    kurtosis /= pow(stdev(input_vec),4.0);
    
    return kurtosis;
}







#endif /* utilities_h */
