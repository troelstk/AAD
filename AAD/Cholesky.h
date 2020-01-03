//
//  Cholesky.h
//  AAD
//
//  Created by Troels Tang Karlsen on 16/12/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef Cholesky_h
#define Cholesky_h

#include <cmath>
  
/* Computes the lower triangular Cholesky decomposition of input matrix A */
template<class T>
vector<vector<T>> Chol( vector<vector<T>> A)
{
    size_t n = A[1].size();
    vector<vector<T>> L(n, vector<T>(n, T(0.0)));
  
    // Decomposing a matrix into Lower Triangular
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            T sum(0.0);
            // Evaluating L(i, j) using L(j, j)
            for(int k = 0; k < j; k++){
                sum += (L[i][k] * L[j][k]);
            }
            L[i][j] = (A[i][j] - sum) / L[j][j];
        }
        T sum(0.0);
        for (int k = 0; k < i; k++) {
            sum += pow(L[i][k], 2.0);
        }
        L[i][i] = sqrt(A[i][i] - sum);
    }
    
    return L;
}

/* Matrix times column vector */
template<class T>
vector<T> simGauss( vector<vector<T>> & L, vector<double> & gauss)
{
    size_t n = gauss.size();
    vector<T> res(n);
    
    for(int i = 0; i<n; i++){
        res[i] = 0.0;
        for(int j = 0; j<i+1; j++){
            res[i] += L[i][j] * gauss[j];
        }
    }
    return res;
}

#endif /* Cholesky_h */
