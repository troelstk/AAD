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
  
template<class T>
vector<vector<T>> myChol( vector<vector<T>> A)
{
    size_t n = A[1].size();
    vector<vector<T>> L(n, vector<T>(n, number(0.0)));
  
    // Decomposing a matrix into Lower Triangular
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            T sum(0.0);
            if (j == i) // summation for diagnols
            {
                for (int k = 0; k < j; k++) {
                    sum += pow(L[j][k], 2.0);
                }
                L[j][j] = sqrt(A[j][j] - sum);
            } else {
                // Evaluating L(i, j) using L(j, j)
                for(int k = 0; k < j; k++){
                    sum += (L[i][k] * L[j][k]);
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    
    // Displaying Lower Triangular and its Transpose
    for (int i = 0; i < n; i++) {
        // Lower Triangular
        for (int j = 0; j < n; j++)
        {
            cout <<  L[i][j].value() << ",";
        }
        cout << "\n";
    }
    
    return L;
}


#endif /* Cholesky_h */
