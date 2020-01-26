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
#include <cassert>
  
/* Computes the lower triangular Cholesky decomposition of input matrix A. Source: https://www.geeksforgeeks.org/cholesky-decomposition-matrix-decomposition/ */
template<class T>
vector<vector<T>> Chol( vector<vector<T>> & A)
{
    size_t n = A[0].size();
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
        // here i==j
        T sum(0.0);
        for (int k = 0; k < i; k++) {
            sum += L[i][k]*L[i][k];
        }
        L[i][i] = sqrt(A[i][i] - sum);
    }
    
    return L;
}

template<class T>
vector<vector<T>> Chol2( vector<vector<T>> & C, vector<vector<T>> & B)
{
    
    int m = int(C.size()); // rows
    int n = int(C[0].size()); // cols
    
    vector<vector<T>> L(m, vector<T>(n, T(0.0)));
    T C_test = T(0.0);
    for(int i=0; i<m; ++i) {

        
        for(int j=i+1; j<n; ++j){
            C_test = 0.0;
            for(int k = 0; k<m; ++k)
            {
                C_test += B[i][k] * B[j][k];
            }
            for(int k = i-1; k >= 0; --k){
                C_test -= L[i][k] * L[j][k];
            }
            if(j==i) {
                L[i][i] = sqrt(C_test);
            }
            else {
                L[j][i] = C_test / L[i][i] ;
            }
        }
        
    }
    return L;
}

template<class T>
void Chol2Adj( vector<vector<T>> & C,
                           vector<vector<T>> & B,
                           vector<vector<T>> & L ) {
    
    int m = C.size(); // rows
    int n = C[0].size(); // cols
    
    for(int i=m-1; i>=0; --i) {
        for(int j=n-1; j>=i; --j){
            if(j==i) {
                C[i][i].adjoint() = 0.5 * L[i][i].adjoint()/L[i][i];
            }
            else {
                C[i][j].adjoint() = L[j][i].adjoint()/L[i][i] ;
                L[i][i].adjoint() -= L[j][i].adjoint()*L[j][i]/L[i][i] ;
            }
            for(int k =0; k<i; ++k){
                L[i][k].adjoint() -= C[i][j].adjoint() * L[j][k];
                L[j][k].adjoint() -= C[i][j].adjoint() * L[i][k];
            }
            for(int k = m-1; k>=0; --k){
                B[i][k].adjoint() -= C[i][j].adjoint() * L[j][k];
                B[j][k].adjoint() -= C[i][j].adjoint() * L[i][k];
            }
        }
    }
}

template<class T> inline
vector<vector<T>> transp(vector<vector<T>> A ) {
    size_t m = A.size(); // rows
    size_t n = A[0].size(); // cols
    
    vector<vector<T>> res(n, vector<T>(m, T(0.0)));
    
    for(size_t i=0; i<m; ++i) {
        for(size_t j=0; j<n; ++j){
            res[i][j] = A[j][i];
        }
    }
    return res;
}

/* Matrix times column vector */
template<class T>
vector<T> simGauss( vector<vector<T>> & L, vector<double> & gauss)
{
    // Check LHS n_cols == RHS n_rows 
    assert( L[0].size() == gauss.size() );
    size_t m = gauss.size();
    vector<T> res(m, T(0.0));
    
    for(int i = 0; i<m; i++){
        for(int j = 0; j<i+1; j++){
            res[i] += L[i][j] * gauss[j];
        }
    }
    return res;
}

/* Matrix times column vector */
template<class T, class S>
vector<T> MatVecProd( vector<vector<T>> & L, vector<S> & gauss)
{
    assert( L[0].size() == gauss.size() );
    size_t n = L.size();
    vector<T> res(n, T(0.0));
    
    for(int i = 0; i<n; i++){
        for(int j = 0; j<L[0].size(); j++){
            res[i] += L[i][j] * gauss[j];
        }
    }
    return res;
}

template<class T>
vector<vector<T>> MatMatProd( vector<vector<T>> LHS, vector<vector<T>> RHS)
{
    size_t LHS_rows = LHS.size();
    size_t LHS_cols = LHS[0].size();  // == RHS_rows
    size_t RHS_rows = RHS.size();
    size_t RHS_cols = RHS[0].size();
    
    assert( LHS_cols == RHS_rows );

    vector<vector<T>> res(LHS_rows, vector<T>(RHS_cols, T(0.0)));
    
    for(int i = 0; i<LHS_rows; i++){
        for(int j = 0; j<RHS_cols; j++){
            for(int k = 0; k<LHS_cols; k++){
                res[i][j] += LHS[i][k] * RHS[k][j];
            }
        }
    }
    return res;
}

#endif /* Cholesky_h */
