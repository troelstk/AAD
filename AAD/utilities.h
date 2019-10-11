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
template<class T>
T sum(vector<T>& input_vec) {
    T res(0.0);
    for( auto &x : input_vec) {
        res += x;
    }
    return res;
}

#endif /* utilities_h */
