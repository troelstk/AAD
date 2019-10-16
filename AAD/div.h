//
//  div.h
//  AAD
//
//  Created by Troels Tang Karlsen on 14/10/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef div_h
#define div_h

// Sample correlation matrix: rho_hat(i,j) = 0.5 + (1-0.5)*exp(-0.05*abs(i-j));
mat rho_hat(10,10);
for(int i = 0; i<10; i++){
    for(int j = 0; j<10; j++){
        rho_hat(i,j) = 0.5 + (1-0.5)*exp(-0.05*abs(i-j));
    }
}
eig_sym(eigval, P, rho_hat);
H = diagmat(eigval);
if( eigval.is_sorted() == 0 ) return 0.0;
P_n = P.tail_cols(n);
largest_n_eigvals = eigval.tail_rows(n);
L_n = sqrt(diagmat(largest_n_eigvals));

B_n = P_n * L_n;
rho_bar_n = B_n*B_n.t();
rho_n = rho_bar_n;

for(int i = 0; i<10; i++){
   for(int j = 0; j<10; j++){
       rho_n(i,j) = rho_bar_n(i,j)/sqrt(rho_bar_n(i,i)*rho_bar_n(j,j));
   }
}
print(rho_n);

#endif /* div_h */
