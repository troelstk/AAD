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
                             int seed1, int nPaths, int nSteps, T yearly_payments, int dim_n)
{
    
    T dt(double(Ta-t)/double(nSteps));
    T tau(1.0/yearly_payments );
    T sqDt = sqrt(dt);
    T res(0.0);
    size_t M1 = initF.size();
    
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

    mat Cs = C( span(int(Ta),M1-1), span(int(Ta),M1-1) ); // Select subset of full covariance matrix to simulate from
    mat corr2 = C2( span(int(Ta),M1-1), span(int(Ta),M1-1) );
    // Find eigenvalues and vectors of correlation matrix
    mat P;
    vec eigval;
    eig_sym(eigval, P, corr2);
    // H has diagonal with eigenvalues
    mat H = diagmat(eigval);
    
    // Choose rank n:
    if( eigval.is_sorted() == 0 ) return 0.0;
    // Assumes eigval is sorted ascending
    vec largest_n_eigvals = eigval.tail_rows(dim_n);
    // H2 is H but only largest n eigenvalues
    mat H_n = diagmat(largest_n_eigvals);
    // Lambda_n = sqrt(H_n)
    mat L_n = sqrt(H_n);
    // P_n is columns corresponding to largest n eigenvalue
    mat P_n = P.tail_cols(dim_n);
    
    // B_n is used to simulate from
    mat B_n = P_n * L_n;
    
    arma_rng::set_seed(seed1);
    
    vector<T> F;
    mrg32k3a myRNG(seed1, 123) ;
    myRNG.init(dim_n);

    vector<double> gauss(dim_n);
    vector<T> lnRates(nPaths);
    
    for(int i = 0; i<nPaths; ++i){
        mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
        vector<T> lnF = log(initF);
        for(int n = 0; n<nSteps; ++n) {
            //myRNG.nextG(gauss);
            //mat gauss2 = B_n * vec(gauss);
            for(int k = int(Ta); k < int(Tb); ++k) // Loop over number of Forward rates to simulate, 0=F1, 9=F10,..,19=F20
            {   // When Ta is 9, we access 10'th (alpha+1) entry in F
                // Compute sum in Brigo 6.53
                T sum(0.0);
                for(int j = int(Ta); j <= k; ++j ) { // Loop over yearly forward rates,
                    // 2.0 = F_alpha+1 = F_2 = k[1], 3.0 = alpha+2 = k[2], 4.0 = alpha+3 = beta = k[3],
                    double Fj = exp(lnF[j]);
                    sum += corr[k][j] * vol[j][j] * Fj / (1.0 + tau * Fj);
                }
                lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_n(k-Ta,n);
                //lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss2(k-Ta);
            }
        }
        // Now have F_k(T_alpha) for k=10,..,M
        F = exp(lnF);
        T floating_swap_rate;
        floating_swap_rate = SR_from_F(F, yearly_payments, (int)(Ta), (int)(Tb) );

        lnRates[i] = log(floating_swap_rate);
        
        T diff_rates;
        diff_rates = max(floating_swap_rate - r_fix, 0.0);  // Payer: floating_swap_rate - r_fix
        
        // Compute value of swap using this floating rate
        double C = C_ab( F, yearly_payments, int(Ta), int(Tb));
        res += C * diff_rates;
    }
    print("std is ", stdev(lnRates));
    print("kurtosis is ", kurtosis(lnRates));
    print("skew is ", skew(lnRates));
    
    // Discounting back to time t:
    T disc;
    disc = DF_from_F(initF, yearly_payments, int(t), int(Ta));
    return disc * res/double(nPaths) * notional;
}


template<class T> T LMM_BermudaSwaption(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF, vector<T> & exTimes,
                             T t, T Ta, T Tb, T r_fix, T notional,
                             int seed1, int nPaths, int nPaths_presim, int nSteps_y, T yearly_payments, int dim_n)
{
    T tau(1.0/yearly_payments );
    
    size_t M1 = initF.size();
    
    vec means(M1, fill::zeros);
    //mat C(M1, M1, fill::zeros);
    mat C2(M1, M1, fill::zeros);
    
    // Fill covariance matrix
    for(int i = 0; i<M1; ++i){
        for(int j = 0; j<M1; ++j){
            //C(i,j) = corr[i][j] * vol[i][i] * vol[j][j];
            C2(i,j) = corr[i][j];
            //print(i, " ", j, " Corr ", corr[i][j], " vol i ", vol[i][i], " vol j ", vol[j][j], " gives ", C(i,j) );
        }
    }

    //mat Cs = C( span(int(Ta),M1-1), span(int(Ta),M1-1) ); // Select subset of full covariance matrix to simulate from
    mat corr2 = C2( span(int(Ta),M1-1), span(int(Ta),M1-1) );
    // Find eigenvalues and vectors of correlation matrix
    mat P;
    vec eigval;
    eig_sym(eigval, P, corr2);
    // H has diagonal with eigenvalues
    mat H = diagmat(eigval);
    
    // Choose rank n:
    if( eigval.is_sorted() == 0 ) return 0.0;
    // Assumes eigval is sorted ascending
    vec largest_n_eigvals = eigval.tail_rows(dim_n);
    // H2 is H but only largest n eigenvalues
    mat H_n = diagmat(largest_n_eigvals);
    // Lambda_n = sqrt(H_n)
    mat L_n = sqrt(H_n);
    // P_n is columns corresponding to largest n eigenvalue
    mat P_n = P.tail_cols(dim_n);
    
    // B_n is used to simulate from
    mat B_n = P_n * L_n;
    
    arma_rng::set_seed(seed1);
    
    vector<T> F;
    mrg32k3a myRNG(seed1, 123) ;
    myRNG.init(dim_n);

    vector<double> gauss(dim_n);
    vector<T> lnRates(nPaths_presim);
    
    vector<vector<double>> SR(nPaths_presim, vector<double>(exTimes.size()));
    vector<vector<double>> Libor(nPaths_presim, vector<double>(exTimes.size()));
    vector<vector<double>> swap_vals(nPaths_presim, vector<double>(exTimes.size()));
    
    for(int t = 1; t<exTimes.size(); ++t ){
        //print("exTime is ", exTimes[t], " t is ", t);
        T swap_val(0.0);
        int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
        T dt( 1.0/double(nSteps));
        T sqDt = sqrt(dt);
        
        for(int i = 0; i<nPaths_presim; ++i){
            mat gauss_n = B_n * arma::mvnrnd(zeros(dim_n), eye(dim_n,dim_n), nSteps);
            vector<T> lnF = log(initF);
            for(int n = 0; n<nSteps; ++n) {
                //print("Time is ", exTimes[t-1] + dt*(n+1));
                for(int k = int(exTimes[t]); k < int(Tb); ++k)
                {
                    T sum(0.0);
                    for(int j = int(exTimes[t]); j <= k; ++j ) {
                        double Fj = exp(lnF[j]);
                        sum += corr[k][j] * vol[j][j] * Fj / (1.0 + tau * Fj);
                    }
                    lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_n(k-Ta,n);
                    //print("F", k+1, " simulated to ", exp(lnF[k]) );
                } // rates
            } // steps
            // Now have F_k(T_alpha) for k=10,..,M
            F = exp(lnF);
            T floating_swap_rate;
            floating_swap_rate = SR_from_F(F, yearly_payments, (int)(exTimes[t]), (int)(Tb) );

            SR[i][t] = floating_swap_rate;
            Libor[i][t] = exp(lnF[exTimes[t]]);
            T diff_rates;
            diff_rates = max(floating_swap_rate - r_fix, 0.0);  // Payer: floating_swap_rate - r_fix
            
            // Compute value of swap using this floating rate
            double C_i = C_ab( F, yearly_payments, int(exTimes[t]), int(Tb));

            swap_val = C_i * diff_rates;
            
            // Save Swap value at time t
            swap_vals[i][t] = swap_val;
        } // paths
        
        // Discounting back to time t:
        T disc;
        disc = DF_from_F(initF, yearly_payments, int(t), int(exTimes[t]));
    } // Exercise dates
    
    // Find regression coefficients
    mat U;
    vec s;
    mat V;
    double lambda = 0.0; // Tikhonov parameter
    int no_basis_funcs = 3;
    
    // Exercise value at last exercise time
    vec Y(nPaths_presim);
    vec payoff(nPaths_presim);
    for(int i = 0; i<nPaths_presim; ++i){
        payoff(i) = swap_vals[i][exTimes.size() - 1];
    }
    
    // Backwards loop:
    for(int t = int(exTimes.size() - 2); t >= 1; --t ){
        //print("exTime is ", exTimes[t], " t is ", t);
        vector<double> ITMswapRate, ITMLibor, ITMY;
        vector<int> indices;
        for(int i = 0; i<nPaths_presim; ++i){
            payoff(i) /=  (1 + Libor[i][t]);
        }
        // Find in the money paths
        for(int i = 0; i<nPaths_presim; ++i){
            if(payoff(i) > 0) { // If ITM
                ITMswapRate.push_back(SR[i][t]);
                ITMLibor.push_back(Libor[i][t]);
                ITMY.push_back(payoff(i));
                indices.push_back(i);
            };
        }
        size_t ITMsize = ITMswapRate.size();
        
        mat X(ITMsize, no_basis_funcs, fill::zeros);

        
        for(int i = 0; i<ITMsize; ++i){
            X(i, 0) = 1;
            X(i, 1) = ITMswapRate[i];
            X(i, 2) = ITMLibor[i];
        }

        // Do SVD:
        svd(U,s,V,X);
        
        mat D(ITMsize, no_basis_funcs, fill::zeros);
        mat Sig(ITMsize, ITMsize, fill::zeros);

        
        for(int i =0; i<s.n_rows; ++i){
            D(i,i) = s(i);
            Sig(i,i) = 1.0/(s(i)*s(i) + lambda*lambda );
        }
        
        vec beta = V * D.t() * Sig * U.t() * vec(ITMY) ;
        
        vec EY = X*beta;
        
        for(int i = 0; i<ITMsize; ++i){
            if( swap_vals[i][t] > EY[i]){
                payoff(indices[i]) = swap_vals[i][t] ;
            }
        }
    } // Exercise times

    double sumPayoff = sum(payoff);
    T disc;
    disc = DF_from_F(initF, yearly_payments, int(t), int(exTimes[1]));
    print("disc is ", disc, " from time ", int(t), " to time ", int(exTimes[1]) );
    return disc * sumPayoff/double(nPaths_presim) * notional;
}


#endif /* LMM_h */

void dummyFunc() {
    // Working american option pricing example from LS:
    int nPaths_presim = 8;
    mat U;
    vec s;
    mat V;
    double lambda = 0.0; // Tikhonov parameter
    int no_basis_funcs = 3;

    // Exercise value at last exercise time

    vec payoff(nPaths_presim);
    mat LSstock = {
      {1.09, 1.08, 1.34},
      {1.16, 1.26, 1.54},
      {1.22, 1.07, 1.03},
      {0.93, 0.97, 0.92},
      {1.11, 1.56, 1.52},
      {0.76, 0.77, 0.90},
      {0.92, 0.84, 1.01},
      {0.88, 1.22, 1.34} };

    print(LSstock);
    for(int i = 0; i<nPaths_presim; ++i){
      payoff(i) = max(1.1 - LSstock(i,2), 0.0);
      print("Payoff at time t= ", 2, " is ", payoff(i));
    }
    // Backwards loop:
    for(int t = 2; t >= 1; --t ){
      vector<double> ITMpaths;
      vector<double> ITMY;
      vector<int> indices;
      for(int i = 0; i<nPaths_presim; ++i){
          payoff(i) /=  (1.06);
      }
      // Find in the money paths
      for(int i = 0; i<nPaths_presim; ++i){
          if(1.1 - LSstock(i,t-1) > 0.0) {
              ITMpaths.push_back(LSstock(i,t-1));
              ITMY.push_back(payoff(i));
              indices.push_back(i);
          };
      }
      size_t ITMsize = ITMpaths.size();
      
      mat X(ITMsize, no_basis_funcs, fill::zeros);
      vec Y(ITMsize);

      
      for(int i = 0; i<ITMsize; ++i){
          X(i, 0) = 1;
          X(i, 1) = ITMpaths[i];
          X(i, 2) = X(i, 1)*X(i, 1);
      }
      // Do SVD:
      svd(U,s,V,X);
      
      mat D(ITMsize, no_basis_funcs, fill::zeros);
      mat Sig(ITMsize, ITMsize, fill::zeros);

      
      for(int i=0; i<s.n_rows; ++i){
          D(i,i) = s(i);
          Sig(i,i) = 1.0/(s(i)*s(i) + lambda*lambda );
      }
      
      vec beta = V * D.t() * Sig * U.t() * vec(ITMY);
      
      vec EY = X * beta;
      
      for(int i = 0; i<ITMsize; ++i){
          if( max(1.1 - ITMpaths[i], 0.0) > EY[i]){
              payoff(indices[i]) = max(1.1 - ITMpaths[i], 0.0) ;
          }
      }

      
    } // Exercise times

    double sumPayoff = arma::sum(payoff)/8.0;

    print("Am option price is ", sumPayoff/1.06);
    
}

/*print("Simulated swap rates");
print(exp(lnRates));
print("mean is ", mean(lnRates));
print("std is ", stdev(lnRates));
print("kurtosis is ", kurtosis(lnRates));
print("skew is ", skew(lnRates));*/
