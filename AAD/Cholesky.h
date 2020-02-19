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
    
    int n = int(C.size()); // rows
    int m = int(B[0].size()); // cols
    
    vector<vector<T>> L(n, vector<T>(m, T(0.0)));
    
    for(int i=0; i<m; ++i) {
        for(int j=i; j<n; ++j){
            C[i][j] = 0.0;
            // Reduced rank step
            for(int k = 0; k<m; ++k)
            {
                C[i][j] += B[i][k] * B[j][k];
            }
            for(int k = i-1; k >= 0; --k){
                C[i][j] -= L[i][k] * L[j][k];
            }
            if(j==i) {
                L[i][i] = sqrt(C[i][j]);
            }
            else {
                L[j][i] = C[i][j] / L[i][i] ;
            }
        }
    }
    return L;
}

template<class T>
void Chol2Adj( vector<vector<T>> & C,
               vector<vector<T>> & B,
               vector<vector<T>> & L ) {
    assert(C.size() == B.size());
    assert(C.size() == L.size());
    assert(B[0].size() == L[0].size());
    int n = int(C.size()); // rows of C == cols of C == rows of B
    int m = int(B[0].size()); // cols of B (dim_n == rank of cov)
    
    //vector<vector<T>> L(n, vector<T>(m, T(0.0)));
    
    for(int i=0; i<n; ++i) {
        for(int j=0; j<m; ++j){
            //L[i][j] = L2[i][j];
        }
    }
    
    //print("n x m ", n, "x", m);
    
    for(int i=m-1; i>=0; --i) {
        for(int j=n-1; j>=i; --j){
            print("pre C[",i, "][",j,"] adj ", C[i][j].adjoint(), " L adj ", L[i][i].adjoint(), " L val ", double(L[i][i]) );
            if(j==i) {
                
                C[i][i].adjoint() = 0.5 * L[i][i].adjoint() / double(L[i][i]);
            }
            else {
                C[i][j].adjoint() = L[j][i].adjoint() / double(L[i][i]);
                L[i][i].adjoint() -= L[j][i].adjoint() * double(L[j][i]) / double(L[i][i]);
            }
            print("post C[",i, "][",j,"] adj ", C[i][j].adjoint(), " L adj ", L[j][i].adjoint(), " L val ", double(L[i][i]) );
            
            for(int k=0; k<i; ++k){
                L[i][k].adjoint() -= C[i][j].adjoint() * double(L[j][k]);
                //print(L[i][k].adjoint());
                L[j][k].adjoint() -= C[i][j].adjoint() * double(L[i][k]);
                //print(L[j][k].adjoint());
            }
            for(int k = m-1; k>=0; --k){
                B[i][k].adjoint() += C[i][j].adjoint() * double(B[j][k]);
                B[j][k].adjoint() += C[i][j].adjoint() * double(B[i][k]);
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
            res[j][i] = A[i][j];
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

/* Matrix times column vector */
template<class T>
vector<vector<T>> MatMatAdd( vector<vector<T>> LHS, vector<vector<T>> RHS )
{
    assert( LHS[0].size() == RHS[0].size() );
    assert( LHS.size() == RHS.size() );
    
    size_t m = LHS.size(); // rows
    size_t n = LHS[0].size(); // cols
    
    vector<vector<T>> res(m, vector<T>(n, T(0.0)));
    
    for(int i = 0; i<m; i++){
        for(int j = 0; j<n; j++){
            res[i][j] = LHS[i][j] + RHS[i][j];
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
