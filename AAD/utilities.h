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
#include "number.h"
static bool debug_var = true;
using namespace std;

// Prints the values of a vector in one row
inline void print(vector<double> input_vec) {
    for( auto &x : input_vec) {
        cout << x << " ";
    }
    cout << endl;
}
template<class number> inline void print(vector<number> input_vec) {
    for( auto &x : input_vec) {
        cout << x.value() << " ";
    }
    cout << endl;
}
// Prints the values of a vector in one row
inline void
print_DEBUG(vector<double> input_vec) {
    if(debug_var){
        print(input_vec);
    }
}
// Prints the values of a vector in one row
template<class number> inline void
print_DEBUG(vector<number> input_vec) {
    if(debug_var){
        print(input_vec);
    }
}

// Prints values of a vector of vectors, in one row per vector
inline void
print(vector<vector<double>> & input_vec_of_vecs) {
    for( auto &x : input_vec_of_vecs) {
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
template<typename... ArgTypes>
void print_DEBUG(ArgTypes... args)
{
    if(debug_var){
        expand_type{ 0, (std::cout << args, 0)... };
        std::cout << endl;
    }
}

// Overload log for vector of T's
template<class T>
vector<T> log(vector<T> & input_vec) {
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

// Mean of a vector
template<class T>
T mean(vector<T>& input_vec) {
    return sum(input_vec)/double(input_vec.size());
}

template<class T>
T stdev(vector<T>& input_vec) {
    T stdev(0.0);
    T Mean = mean(input_vec);
    
    for(auto &x : input_vec){
        stdev += pow(x - Mean, 2.0) ;
    }
    stdev /= double(input_vec.size());
    
    return sqrt(stdev);
}

template<class T>
T skew(vector<T>& input_vec) {
    T skew(0.0), skewtop(0.0), skewbot(0.0);
    T Mean = mean(input_vec);
    
    for(auto &x : input_vec){
        skewtop +=  pow(x - Mean, 3.0)/double(input_vec.size());
        skewbot +=  pow(x - Mean, 2.0)/(double(input_vec.size())-1.0);
    }
    skew = skewtop/ pow(skewbot, 1.5);
    
    return skew;
}

template<class T>
T kurtosis(vector<T>& input_vec) {
    T kurtosis(0.0);
    T Mean = mean(input_vec);
    
    for(auto &x : input_vec){
        kurtosis += pow(x - Mean, 4.0)/double(input_vec.size());
    }
    kurtosis /= pow(stdev(input_vec),4.0);
    
    return kurtosis;
}


inline void writeToFile(string MyfileName, arma::mat inpMat ){
    ofstream myfile (MyfileName);
    if (myfile.is_open())
    {
        for(int i = 0; i<inpMat.n_rows; ++i){
            for(int j = 0; j<inpMat.n_cols; ++j ){
                myfile << inpMat(i, j) << ",";
            }
            myfile << "\n";
        }
        myfile.close();
    }
    else cout << "Unable to open file\n";
}



#endif /* utilities_h */
