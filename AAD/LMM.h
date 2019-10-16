//
//  LMM.h
//  AAD
//
//  Created by Troels Tang Karlsen on 22/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef LMM_h
#define LMM_h
#include <cmath>
#include "utilities.h"
#include "gaussians.h"
#include "mrg32k3a.h"
#include "/usr/local/Cellar/armadillo/9.700.2/include/armadillo"

using namespace std;

template<class T> T LMM_swaption(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF,
                             T t, T Ta, T Tb, T r_fix, T notional,
                             int seed1, int nPaths, int nSteps, T yearly_payments)
{
    
    T dt(double(Ta-t)/double(nSteps));
    T tau(1.0/yearly_payments );
    T sqDt = sqrt(dt);
    T res(0.0);
    T avg_float(0.0);
    size_t M1 = initF.size();
    
    //T nPayments( (Tb - Ta) * yearly_payments);
    //print("payments ", nPayments, " ", M);
    
    vec means(M1, fill::zeros);
    mat C(M1, M1, fill::zeros);
    mat C2(M1, M1, fill::zeros);
    
    // Fill covariance matrix
    for(int i = 0; i<M1; ++i){
        for(int j = 0; j<M1; ++j){
            C(i,j) = corr[i][j] * vol[i][i] * vol[j][j];
            C2(i,j) = corr[i][j];
            //print(i, " ", j, " Corr ", corr[i][j], " vol i ", vol[i][i], " vol j ", vol[j][j], " gives ", C(i,j) );
        }
    }
    //print("Cov is"); print(C);
    mat Cs = C( span(int(Ta),M1-1), span(int(Ta),M1-1) ); // Select subset of full covariance matrix to simulate from
    mat corr2 = C2( span(int(Ta),M1-1), span(int(Ta),M1-1) );
    // Find eigenvalues and vectors of correlation matrix
    mat P;
    vec eigval;
    eig_sym(eigval, P, corr2);
    // H has diagonal with eigenvalues
    mat H = diagmat(eigval);
    
    // Choose rank n:
    int n = 2;
    if( eigval.is_sorted() == 0 ) return 0.0;
    // Assumes eigval is sorted ascending
    vec largest_n_eigvals = eigval.tail_rows(n);
    // H2 is H but only largest n eigenvalues
    mat H_n = diagmat(largest_n_eigvals);
    // Lambda_n = sqrt(H_n)
    mat L_n = sqrt(H_n);
    // P_n is columns corresponding to largest n eigenvalue
    mat P_n = P.tail_cols(n);
    
    // B_n is used to simulate from
    mat B_n = P_n * L_n;
    print(size(B_n));
    print(B_n);
    
    
    
    //print(size(Cs));
    //print("Cov small is"); print(Cs);
    //vec meansS(M1, fill::zeros);
    //vec means(n, fill::zeros);
    
    arma_rng::set_seed(seed1);
    
    vector<T> F;
    int count = 0;
    mrg32k3a myRNG(seed1, 123) ;
    myRNG.init(n);

    vector<double> gauss(n);
    
        
    vector<T> lnRates(nPaths);
    //mat myMat = mvn_rnd(means, Cs, 1);
    
    for(int i = 0; i<nPaths; ++i){
        // Simulate multivariate gaussians with covariance C and mean 0
        //mat gauss = arma::mvnrnd(means, C, nSteps);
        //mat gaussS = arma::mvnrnd(meansS, Cs, nSteps);
        mat gauss_n = B_n * arma::mvnrnd(zeros(n), eye(n,n), nSteps);
        
        
        vector<T> lnF = log(initF);
        for(int n = 0; n<nSteps; ++n) {
           // myRNG.nextG(gauss);
           // vec arma_gauss = gauss;
           // mat gauss_new = B_n * arma_gauss;
            //print("STEP ", n);
            for(int k = int(Ta); k < int(Tb); ++k) // Loop over number of Forward rates to simulate, 0=F1, 9=F10,..,19=F20
            {   // When Ta is 9, we access 10'th (alpha+1) entry in F
                //print("k is ", k);
                // Compute sum in Brigo 6.53
                T sum(0.0);
                for(int j = int(Ta); j <= k; ++j ) { // Loop over yearly forward rates,
                    // 2.0 = alpha+1 = k[1], 3.0 = alpha+2 = k[2], 4.0 = alpha+3 = beta = k[3],
                    //print("j is ", j);
                    double Fj = exp(lnF[j]);
                    //print(Fj);
                    sum += corr[k][j] * vol[j][j] * Fj / (1.0 + tau * Fj);
                }
                //lnF[k] += vol[k][k] * tau * sum * dt - vol[k][k] * vol[k][k]/2.0 * dt + sqDt * gaussS(k-Ta,n);
                lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_n(k-Ta,n);
                //lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_new(k-Ta);
            }
            //print(exp(lnF));
        }
        // Now have F_k(T_alpha) for k=10,..,M
        F = exp(lnF);
        //print("Sim F is "); print(F);
        T floating_swap_rate;
        floating_swap_rate = SR_from_F(F, yearly_payments, (int)(Ta), (int)(Tb) );
        //print("floating_swap_rate is ", floating_swap_rate);
        //print(floating_swap_rate);
        lnRates[i] = log(floating_swap_rate);
        
        T diff_rates;
        diff_rates = max(floating_swap_rate - r_fix, 0.0);  // Payer: floating_swap_rate - r_fix
        //print("diff_rates is ", diff_rates);
        count += diff_rates > 0.0 ? 1 : 0;
        
        // Compute value of swap using this floating rate
        T val(0.0);
        for(int j = int(Ta) + 1; j <= Tb; ++j){ // Time 10... 20 back to time 9
            // Discount back to time Ta
            T disc = DF_from_F(F, yearly_payments, int(Ta), j);
            //print("Disc is ", disc, " from time ", int(Ta), " to time ",  j);
            val += disc;
            //print(val);
        }
        res += val * diff_rates;
        //print(val * diff_rates * notional);
        avg_float += floating_swap_rate;
    }
    double avg = avg_float/double(nPaths);
    print("avg_float is ", avg);
    double lnMean = sum(lnRates)/double(nPaths);
    
    print("count is ", count);
    
    double std = 0;
    double kurtosis = 0;
    double skewtop = 0.0;
    double skewbot = 0.0;
    
    for(int i = 0; i<nPaths; i++){
        std +=      pow(lnRates[i] - lnMean, 2.0)/nPaths;
        kurtosis += pow(lnRates[i] - lnMean, 4.0)/nPaths;
        skewtop +=  pow(lnRates[i] - lnMean, 3.0)/nPaths;
        skewbot +=  pow(lnRates[i] - lnMean, 2.0)/(nPaths-1);
    }
    double skew = skewtop/ pow(skewbot, 1.5);
    std = sqrt(std);
    kurtosis /= pow(std,4.0);

    print("std is ", std);
    print("kurtosis is ", kurtosis);
    print("skew is ", skew);

    //print(rates);
    
    // Discounting back to time t:
    T disc;
    disc = DF_from_F(initF, yearly_payments, int(t), int(Ta));
    //print("Disc is ", disc);
    return disc * res/double(nPaths) * notional;
}


#endif /* LMM_h */


/*double testprod = 1.0, testprod2 = 1.0, sum = 0.0;
for(int h = int(Ta); h<Tb; ++h) {
    testprod *= 1.0/(1.0+tau*F[h]);
}
for(int k = int(Ta) ; k<Tb; ++k){
    testprod2 = 1.0;
    for(int h = int(Ta); h<=k; ++h) {
        testprod2 *= 1.0/(1.0+tau*F[h]);
    }
    sum += testprod2;
}
double swap_r = (1 - testprod)/sum;
print("new float ", swap_r);*/
/*double prod = 1.0;
for(int i = int(Ta); i<j; ++i){
    prod *= 1.0/(1.0 + tau*F[i]);
}
print("prod is ", prod);*/

// New rank-n approximated correlation matrix rho_n
/*mat rho_bar_n = B_n*B_n.t();
mat rho_n = rho_bar_n;
int M = Tb-Ta;
for(int i = 0; i<M; i++){
    for(int j = 0; j<M; j++){
        rho_n(i,j) = rho_bar_n(i,j)/sqrt(rho_bar_n(i,i)*rho_bar_n(j,j));
    }
}*/
//print(rho_n);
