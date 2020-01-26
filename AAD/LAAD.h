//
//  LAAD.h - Linear Algebra Adjoint Differentiation
//  AAD
//
//  Created by Troels Tang Karlsen on 25/01/2020.
//  Copyright © 2020 Troels Tang Karlsen. All rights reserved.
//

#ifndef LAAD_h
#define LAAD_h

using namespace std

template <typename T>
class Mat {

private:
    size_t n_cols;
    size_t n_rows;
    
    vector<vector<T>> data;

public:
    Mat(int m = 0, int n = 0 ) : m(n_rows), n(n_cols)  {
        vector<vector<T>> new_data(m, vector<T>(n));
        data = new_data;
    }
    

    // -
    
    // *
    
    // +
    
    // hadamard
    
    // .t()
    
    // .tail_cols(n)
    
    // .swap_cols()
    
    explicit operator T& (int i, int j)
    {
        return data[i][j];
    };
    
}


#endif /* LAAD_h */
