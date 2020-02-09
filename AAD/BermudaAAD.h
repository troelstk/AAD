//
//  BermudaAAD.h
//  AAD
//
//  Created by Troels Tang Karlsen on 18/11/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef BermudaAAD_h
#define BermudaAAD_h

#include <cmath>
#include "utilities.h"
#include "gaussians.h"
#include "mrg32k3a.h"
#include "/usr/local/Cellar/armadillo/9.700.2/include/armadillo"
#include "Cholesky.h"

template<class T> T LMM_BermudaSwaptionAAD(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF, vector<double> & exTimes,
                             double t, double Ta, double Tb, T r_fix, double notional,
                             int seed1, int seed2, int nPaths, int nPaths_presim, int nSteps_y, double yearly_payments, int dim_n)
{
    T tau(1.0/yearly_payments );
    int int_Ta = (int)(Ta), int_Tb = (int)(Tb);
    int no_basis_funcs = 3;
    mat beta(no_basis_funcs, exTimes.size() - 1);
    double lambda = 0.0; // Tikhonov parameter
    
    size_t M1 = initF.size();
    size_t M = int_Tb-int_Ta;
    
    //double eps = 10e-10; // Added to diagonal to make matrix positive definite
    
    /*vector<vector<T>> cov_s(M, vector<T>(M));
    vector<vector<T>> lower(M, vector<T>(M));
    vector<vector<double>> cov_d(  M, vector<double>(M));
    vector<vector<double>> lower_d(M, vector<double>(M));
    
    // Fill covariance vec of vecs
    for(int i = int_Ta; i<Tb; ++i){
        for(int j = int_Ta; j<Tb; ++j){
            T cov_ij(corr[i][j] * vol[i][i] * vol[j][j]);
            cov_s[i-int_Ta][j-int_Ta] = cov_ij ;
            cov_d[i-int_Ta][j-int_Ta] = cov_ij.value();
        }
        cov_s[i-int_Ta][i-int_Ta] += eps;
        cov_d[i-int_Ta][i-int_Ta] += eps;
    }
    
   
    lower = Chol(cov_s);
    lower_d = Chol(cov_d);*/
    
    vector<vector<T>> Cov_num(M, vector<T>(M));
    mat Cov(M, M, fill::zeros);
    
    // Fill covariance matrix
    for(int i = int_Ta; i<int_Tb; ++i){
        int i2 = i-int_Ta;
        for(int j = int_Ta; j<int_Tb; ++j){
            
            int j2 = j-int_Ta;
            Cov_num[i2][j2] = corr[i][j]; //*vol[i][i]*vol[j][j];
            Cov(i2,j2) = double(Cov_num[i2][j2]);
        }
        //Cov(i2,i2) += eps;
    }
    
    mat P; // U
    vec eigval; // D
    eig_sym(eigval, P, Cov);
    //print(eigval);
    //print(P);
    vec largest_n_eigvals = eigval.tail_rows(dim_n);
    mat H_n = diagmat(largest_n_eigvals);
    mat L_n = sqrt(H_n);
    mat P_n = P.tail_cols(dim_n);
    mat B_n = P_n * L_n;
    
    vector<vector<T>> L_num(dim_n, vector<T>(dim_n, T(0.0)));
    vector<vector<T>> eig_num(dim_n, vector<T>(dim_n, T(0.0)));
    vector<vector<T>> P_num(M, vector<T>(dim_n, T(0.0)));
    vector<vector<T>> B_num(M, vector<T>(dim_n, T(0.0)));
    
    for(int i = 0; i<dim_n; ++i){
        eig_num[i][i] = T(largest_n_eigvals(i));
        L_num[i][i] = sqrt(eig_num[i][i]);
    }
    
    // Copy eigen vectors of largest eigenvalues to overloaded P
    for(int i = 0; i<M; i++){ // Rows
        for(int j = 0; j<dim_n; j++){ // cols
            P_num[i][j] = T(P_n(i,j));
        }
    }

    B_num = MatMatProd(P_num, L_num);
    
    
    // Pre compute product, instead of doing it in inner most loop
    vector<vector<T>> corr_vol(M1, vector<T>(M1));
    for(int k = 1; k < int_Tb; ++k)
    {
        for(int j = 1; j <= k; ++j ) {
            corr_vol[k][j] = corr[k][j] * vol[j][j];
        }
    }
    
    vector<double> Fdouble;

    clock_t begin_time = clock();
    
    T disc;
    int int_t = int(t), int_ExT1 = int(exTimes[1]);
    disc = DF_from_F(initF, yearly_payments, int_t, int_ExT1);
    
    vector<double> initFdouble;
    for(auto & x : initF ) {initFdouble.push_back(double(x));}
    
    
    if(exTimes.size() > 2)
    { // Scope to destruct everything declared in pre-simulation and backward loop
        mrg32k3a myRNG(seed1, seed2);
        myRNG.init(dim_n);
        vector<double> gauss(dim_n);
        vector<vector<double>> SR(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> Libor(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> swap_vals(nPaths_presim, vector<double>(exTimes.size()));
        
        // Simulate swap rates, pre-simulation
        double dt( 1.0 /double(nSteps_y));
        double sqDt = sqrt(dt);

        for(int i = 0; i<nPaths_presim; ++i) {
            for(int t = 1; t<exTimes.size(); ++t ) {
                //print("exTime is ", exTimes[t], " t is ", t);
                double swap_val(0.0);
                int int_ExTime = (int)(exTimes[t]);
                int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
                vector<double> lnF = log(initFdouble);
                
                for(int n = 0; n<nSteps; ++n) {
                    myRNG.nextG(gauss);
                    //vec gauss_n = B_n * vec(gauss);
                    vec gauss_corr = B_n * vec(gauss);
                    //vector<double> gauss_cov = simGauss( lower_d, gauss);
                    //print(gauss_n);
                    //print_DEBUG("Time is ", exTimes[t-1] + dt*(n+1));
                    for(int k = int_ExTime; k < int_Tb; ++k)
                    {
                        double sum(0.0);
                        for(int j = int_ExTime; j <= k; ++j ) {
                            double Fj = exp(lnF[j]);
                            sum += double(corr_vol[k][j]) * Fj / (1.0 + double(tau) * Fj);
                        }
                        // For cholesky of Cov
                        //lnF[k] += vol[k][k].value()*tau.value()*sum*dt - vol[k][k].value()*vol[k][k].value()/2.0*dt + sqDt*gauss_cov[k-int_Ta];
                        // For PCA of Cov
                        //lnF[k] += vol[k][k].value()*tau.value()*sum*dt - vol[k][k].value()*vol[k][k].value()/2.0*dt + sqDt*gauss_n[k-int_Ta];
                        // For PCA of Corr
                        lnF[k] += double(vol[k][k])*double(tau)*sum*dt - double(vol[k][k])*double(vol[k][k])/2.0*dt + sqDt*double(vol[k][k])*gauss_corr[k-int_Ta];
                    } // rates
                } // steps
                // Now have F_k(T_alpha) for k=10,..,M
                Fdouble = exp(lnF);
                double floating_swap_rate = SR_from_F(Fdouble, yearly_payments, int_ExTime, int_Tb );
                
                swap_val = double(disc) * notional *
                    C_ab( Fdouble, yearly_payments, int_ExTime, int_ExTime, int_Tb) *
                    max(floating_swap_rate - double(r_fix), 0.0);
                SR[i][t] = floating_swap_rate;
                Libor[i][t] = Fdouble[ int_ExTime ];
                swap_vals[i][t] = swap_val;
                
            } // exTimes
        } // paths
        print_DEBUG("Presim forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
        begin_time = clock();
        
        // Find regression coefficients in backwards loop
        // Exercise value at last exercise time
        vec Y(nPaths_presim);
        vec payoff(nPaths_presim);
        for(int i = 0; i<nPaths_presim; ++i){
            payoff(i) = swap_vals[i][exTimes.size() - 1];
        }
        
        // Backwards loop:
        for(int t = int(exTimes.size() - 2); t >= 1; --t )
        {
            //print("exTime is ", exTimes[t], " t is ", t);
            vector<double> ITMswapRate, ITMLibor, ITMY;
            vector<int> indices;
            for(int i = 0; i<nPaths_presim; ++i){
                payoff(i) /=  (1.0 + Libor[i][t]);
            }
            // Find in the money paths
            for(int i = 0; i<nPaths_presim; ++i){
                if(swap_vals[i][t] > 0) { // If ITM
                    ITMswapRate.push_back(SR[i][t]);
                    ITMLibor.push_back(Libor[i][t]);
                    ITMY.push_back(payoff(i));
                    indices.push_back(i);
                };
            }
            size_t ITMsize = ITMswapRate.size();
            mat X(ITMsize, no_basis_funcs, fill::zeros);

            for(int i = 0; i<ITMsize; ++i){
                X.row(i) = vec({ 1.0, ITMswapRate[i], ITMLibor[i] }).t();
            }
            // Do SVD:
            mat U, V;
            vec s;
            svd(U,s,V,X);
            
            mat D(ITMsize, no_basis_funcs, fill::zeros);
            mat Sig(ITMsize, ITMsize, fill::zeros);

            for(int i=0; i<s.n_rows; ++i){
                D(i,i) = s(i);
                Sig(i,i) = 1.0/(s(i)*s(i) + lambda*lambda );
            }
            
            beta.col(t) = V * D.t() * Sig * U.t() * vec(ITMY) ;
        } // Backward loop
        print_DEBUG("Presim backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
        begin_time = clock();
    } // Scope
    
    
    
    // Main-simulation
    vector<T> F;

    T dt( 1.0/double(nSteps_y) );
    T sqDt = sqrt(dt);

    // Main-simulation
    T end_result(0.0);
    
    vector<double> swap_rates(nPaths);
    mrg32k3a myRNG(seed2, seed1);
    myRNG.init(dim_n);
    vector<double> gauss(dim_n);
    
    number::tape->mark();
    for(int i = 0; i<nPaths; ++i){
        number::tape->rewindToMark();
        double eta = 1;
        T swap_path(0.0);
        for(int t = 1; t<exTimes.size(); ++t ){
            //print("exTime is ", exTimes[t], " t is ", t);
            int int_ExTime = (int)(exTimes[t]);
            int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
            T swap_val(0.0);
            vector<T> lnF = log(initF);
            for(int n = 0; n<nSteps; ++n) {
                myRNG.nextG(gauss);
                //vector<T> gauss_n = MatVecProd( B_num, gauss);
                vector<T> gauss_corr = MatVecProd( B_num, gauss);
                for(int k = int_ExTime; k < int_Tb; ++k)
                {
                    T sum(0.0);
                    for(int j = int_ExTime; j <= k; ++j ) {
                        T Fj(exp(lnF[j]));
                        sum += corr_vol[k][j] * Fj / (1.0 + tau * Fj);
                    }
                    //lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*gauss_n[k-int_Ta];
                    lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_corr[k-int_Ta];
                    //print(lnF[k].value());
                } // rates
            } // steps
            // Now have F_k(T_alpha) for k=10,..,M
            F = exp(lnF);

            T floating_swap_rate;
            floating_swap_rate = SR_from_F(F, yearly_payments, int_ExTime, int_Tb );

            swap_val = disc * notional *
                C_ab( F, yearly_payments, int_ExTime, int_ExTime, int_Tb) *
                max(floating_swap_rate - r_fix, 0.0);

            swap_path += eta * swap_val;
            
            double EY2 = t < exTimes.size() - 1 ?
                // continuation value - at last ex time the eta value computed below will never be used (so 0.0 could be anything)
                as_scalar(vec( { 1.0, double(floating_swap_rate), double(F[ int_ExTime ]) }).t() * beta.col(t)) : 0.0;
            eta *= EY2 > swap_val ? 1.0 : 0.0; // 1 if continue, 0 if exercise
            //print("Swap value ", swap_val.value());
        }
        //print("Swap path is ", swap_path.value());
        swap_path.propagateToMark();
        //print(swap_path.value());
        end_result += swap_path;
    }
    print_DEBUG("Main forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    
    number::propagateMarkToStart();
    
    // Propagate B_num to its children
    
    mat D_bar(P.n_rows, P.n_cols, fill::zeros);
    mat P_bar(P.n_rows, P.n_cols, fill::zeros);
    mat F_mat(P.n_rows, P.n_cols, fill::zeros);
    
    // Fill last dim_n columns of P_bar (= U_bar)
    for(size_t i = 0; i<P.n_rows; i++ ){
        for(size_t j = M-dim_n; j < M; ++j ){
            size_t j2 = j-(M-dim_n) ;
            P_bar(i,j) = P_num[i][j2].adjoint();
        }
    }
    
    // Fill adjoints of eig_num = D_bar
    for(size_t i = M-dim_n; i < M; ++i ){
        size_t i2 = i-(M-dim_n);
        D_bar(i,i) = eig_num[i2][i2].adjoint();
    }

    
    // Fill F matrix
    for(size_t i = 0; i<dim_n; ++i ){
        for(size_t j = 0; j < dim_n; ++j ){
            if( i != j) {
                //F_mat(i,j) = 1.0/(D_bar(j,j)-D_bar(i,i));
                F_mat(i,j) = 1.0/(largest_n_eigvals(j)-largest_n_eigvals(i));
            }
        }
    }
    
     // U == P
    mat A_bar = P.i().t()*(D_bar + F_mat % P.t() * P_bar ) * P.t();
    
    /*
    print("P_bar");
    print(P_bar);
    print("D_bar");
    print(D_bar);
    print("F_mat");
    print(F_mat);*/
    print("A_bar");
    print(A_bar);
    
    for(int i = int_Ta; i<int_Tb; ++i){
        for(int j = int_Ta; j<int_Tb; ++j){
            int i2 = i-int_Ta;
            int j2 = j-int_Ta;
            //print("Cov num pre ", Cov_num[i2][j2].adjoint() );
            Cov_num[i2][j2].adjoint() += A_bar(i2,j2);
            //print("Cov num post ", Cov_num[i2][j2].adjoint() );
            //print("i2 ", i2, " j2 ", j2);
            //print("Cov Adj ", vol[i][j].adjoint());
            //print("Corr Adj ", corr[i][j].adjoint());
            
            // Propagate Cov_num node
            Cov_num[i2][j2].getNode()->propagateOne();
            
            // Propagate previous node
            //prev(Cov_num[i2][j2].getNode())->propagateOne();
            
            /*auto prevNode = Cov_num[i2][j2].getNode();
            for(int k = 0; k<1; k++){
                prevNode = prev(prevNode);
                prevNode->propagateOne();
            }*/
            //print("Cov Adj2 ", vol[i][j].adjoint());
            //print("Corr Adj2 ", corr[i][j].adjoint());
        }
    }
    


    return end_result/double(nPaths) ;
}






template<class T> T LMM_BermudaSwaptionAADChol(vector<vector<T>> & vol, vector<vector<T>> & corr,
                             vector<T> & initF, vector<double> & exTimes,
                             double t, double Ta, double Tb, T r_fix, double notional,
                             int seed1, int seed2, int nPaths, int nPaths_presim, int nSteps_y, double yearly_payments, int dim_n)
{
    T tau(1.0/yearly_payments );
    int int_Ta = (int)(Ta), int_Tb = (int)(Tb);
    int no_basis_funcs = 3;
    mat beta(no_basis_funcs, exTimes.size() - 1);
    double lambda = 0.0; // Tikhonov parameter
    
    size_t M1 = initF.size();
    size_t M = int_Tb-int_Ta;
    
    vector<vector<T>> cov_s(M, vector<T>(M));
    vector<vector<T>> lower(M, vector<T>(dim_n));
    vector<vector<double>> cov_d(M, vector<double>(M));
    vector<vector<double>> lower_d(M, vector<double>(dim_n));
    
    mat Cov(M, M, fill::zeros);
    // Fill covariance vec of vecs
    for(int i = int_Ta; i<Tb; ++i){
        for(int j = int_Ta; j<Tb; ++j){
            //cov_s[i-int_Ta][j-int_Ta] = corr[i][j] * vol[i][i] * vol[j][j] ;
            //T cov_ij(corr[i][j] * vol[i][i] * vol[j][j]);
            
            T cov_ij(corr[i][j]);
            cov_s[i-int_Ta][j-int_Ta] = cov_ij;
            
            cov_d[i-int_Ta][j-int_Ta] = double(cov_ij);
            Cov(i-int_Ta, j-int_Ta) = double(cov_ij);
        }
    }
    mat P; // U
    vec eigval; // D
    eig_sym(eigval, P, Cov);
    //print(eigval);
    //print(P);
    vec largest_n_eigvals = eigval.tail_rows(dim_n);
    mat H_n = diagmat(largest_n_eigvals);
    mat L_n = sqrt(H_n);
    mat P_n = P.tail_cols(dim_n);
    mat B_n = P_n * L_n;
    /*
    mat approxCov = B_n * B_n.t();
    print("B*Bt is ");
    print(approxCov);*/
    
    vector<vector<T>> Pert =
      {{T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) },
       {T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) },
       {T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0) },
       {T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) },
       {T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0) },
       {T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) },
       {T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(0.0) },
       {T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0) },
       {T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) },
       {T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0) },
       {T(0.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) }};
    
    vector<vector<T>> L_num(dim_n, vector<T>(dim_n, T(0.0)));
    vector<vector<T>> eig_num(dim_n, vector<T>(dim_n, T(0.0)));
    vector<vector<T>> P_num(M, vector<T>(dim_n, T(0.0))); // (rows, vector<T>(cols, T(0.0)));
    vector<vector<T>> B_num(M, vector<T>(dim_n, T(0.0)));

    
    for(int i = 0; i<dim_n; ++i){
        eig_num[i][i] = T(largest_n_eigvals(i));
        L_num[i][i] = sqrt(eig_num[i][i]);
    }
    
    // Copy eigen vectors of largest eigenvalues to overloaded P
    for(int i = 0; i<M; i++){ // Rows
        for(int j = 0; j<dim_n; j++){ // cols
            P_num[i][j] = T(P_n(i,j));
        }
    }

    B_num = MatMatProd(P_num, L_num);
    
    //vector<vector<T>> test = MatMatProd(MatMatProd(Pert, cov_s), transp(Pert));
    //print("Test P cov P^t");
    //print_val(test);
    
    /*vector<vector<T>> lower2(M, vector<T>(M));
    lower2 = Chol2(cov_s, B_num);

    print("Chol2 of cov_s");
    print_val(lower2);
    
    vector<vector<T>> BBT = MatMatProd(B_num, transp(B_num));
    print("B*B^T is ");
    print_val(BBT); // == cov_s
    print("cov_s");
    print_val(cov_s);*/
    
    /*lower = Chol2(cov_s, B_num);
    print("Chol2 of B*B^T");
    print_val(lower);
    
    vector<vector<T>> BBT2 = MatMatProd(lower, transp(lower));
    print("Test if L*L^T from BBT yields cov_s ");
    print_val(BBT2);
    */
    vector<vector<T>> A = MatMatProd(MatMatProd(transp(Pert), cov_s), Pert);
    print("A= P^t*C*P");
    print_val(A);
    //print("eig vecs");
    //print(P_n);

    lower = Chol2(A, B_num);

    print("Chol of A is ");
    print_val(lower);
    
    /*vector<vector<T>> B_num_swap(M, vector<T>(dim_n, T(0.0)));
    for(int i=0; i<M; ++i) {
        B_num_swap[i][1] = B_num[i][0];
        B_num_swap[i][0] = B_num[i][1];
    }*/
    
    for(int i = 0; i<M; i++){ // Rows
        for(int j = 0; j<dim_n; j++){ // cols
            lower_d[i][j] = double(lower[i][j]);
        }
    }
    
    vector<vector<T>> corr_vol(M1, vector<T>(M1));
    // Pre compute product, instead of doing it in inner most loop
    for(int k = 1; k < int_Tb; ++k)
    {
        for(int j = 1; j <= k; ++j ) {
            corr_vol[k][j] = corr[k][j] * vol[j][j];
        }
    }
    
    vector<double> Fdouble;

    clock_t begin_time = clock();
    
    T disc;
    int int_t = int(t), int_ExT1 = int(exTimes[1]);
    disc = DF_from_F(initF, yearly_payments, int_t, int_ExT1);
    
    vector<double> initFdouble;
    for(auto & x : initF ) {initFdouble.push_back(double(x));}
    
    { // Scope to destruct everything declared in pre-simulation and backward loop
        mrg32k3a myRNG(seed1, seed2);
        myRNG.init(dim_n);
        vector<double> gauss(dim_n);
        //myRNG.init(M);
        //vector<double> gauss(M);
        vector<vector<double>> SR(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> Libor(nPaths_presim, vector<double>(exTimes.size()));
        vector<vector<double>> swap_vals(nPaths_presim, vector<double>(exTimes.size()));
        
        // Simulate swap rates, pre-simulation
        double dt( 1.0 /double(nSteps_y));
        double sqDt = sqrt(dt);

        for(int i = 0; i<nPaths_presim; ++i) {
            for(int t = 1; t<exTimes.size(); ++t ) {
                //print("exTime is ", exTimes[t], " t is ", t);
                
                // Inputs: t, exTimes, nSteps_y, initFdouble, lower, corr_vol, tau, vol, r_fix, int_Tb, notional, disc
                // Outputs: swap_val, floating_swap_rate, Fdouble
                double swap_val(0.0);
                int int_ExTime = (int)(exTimes[t]);
                int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
                vector<double> lnF = log(initFdouble);
                
                for(int n = 0; n<nSteps; ++n) {
                    myRNG.nextG(gauss);
                    //vec gauss_n = B_n * vec(gauss);
                    vec gauss_corr = B_n * vec(gauss);
                    //vector<double> gauss_cov = MatVecProd( lower_d, gauss);
                    //print(gauss_cov);
                    //print_DEBUG("Time is ", exTimes[t-1] + dt*(n+1));
                    for(int k = int_ExTime; k < int_Tb; ++k)
                    {
                        double sum(0.0);
                        for(int j = int_ExTime; j <= k; ++j ) {
                            double Fj = exp(lnF[j]);
                            sum += double(corr_vol[k][j]) * Fj / (1.0 + double(tau) * Fj);
                        }
                        // For cholesky and PCA of Cov
                        //lnF[k] += double(vol[k][k])*double(tau)*sum*dt - double(vol[k][k])*double(vol[k][k])/2.0*dt + sqDt*gauss_cov[k-int_Ta];
                        // For PCA of Cov
                        //lnF[k] += double(vol[k][k])*double(tau)*sum*dt - double(vol[k][k])*double(vol[k][k])/2.0*dt + sqDt*gauss_n[k-int_Ta];
                        // For PCA of Corr
                        lnF[k] += double(vol[k][k])*double(tau)*sum*dt - double(vol[k][k])*double(vol[k][k])/2.0*dt + sqDt*double(vol[k][k])*gauss_corr[k-int_Ta];
                    } // rates
                } // steps
                // Now have F_k(T_alpha) for k=10,..,M
                Fdouble = exp(lnF);
                double floating_swap_rate = SR_from_F(Fdouble, yearly_payments, int_ExTime, int_Tb );
                
                swap_val = double(disc) * notional *
                    C_ab( Fdouble, yearly_payments, int_ExTime, int_ExTime, int_Tb) *
                    max(floating_swap_rate - double(r_fix), 0.0);
                SR[i][t] = floating_swap_rate;
                Libor[i][t] = Fdouble[ int_ExTime ];
                swap_vals[i][t] = swap_val;
            } // exTimes
        } // paths
        print_DEBUG("Presim forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
        begin_time = clock();
        
        // Find regression coefficients in backwards loop
        // Exercise value at last exercise time
        vec Y(nPaths_presim);
        vec payoff(nPaths_presim);
        for(int i = 0; i<nPaths_presim; ++i){
            payoff(i) = swap_vals[i][exTimes.size() - 1];
        }
        
        // Backwards loop:
        for(int t = int(exTimes.size() - 2); t >= 1; --t )
        {
            //print("exTime is ", exTimes[t], " t is ", t);
            vector<double> ITMswapRate, ITMLibor, ITMY;
            vector<int> indices;
            for(int i = 0; i<nPaths_presim; ++i){
                payoff(i) /=  (1.0 + Libor[i][t]);
            }
            // Find in the money paths
            for(int i = 0; i<nPaths_presim; ++i){
                if(swap_vals[i][t] > 0) { // If ITM
                    ITMswapRate.push_back(SR[i][t]);
                    ITMLibor.push_back(Libor[i][t]);
                    ITMY.push_back(payoff(i));
                    indices.push_back(i);
                };
            }
            size_t ITMsize = ITMswapRate.size();
            mat X(ITMsize, no_basis_funcs, fill::zeros);

            for(int i = 0; i<ITMsize; ++i){
                X.row(i) = vec({ 1.0, ITMswapRate[i], ITMLibor[i] }).t();
            }
            // Do SVD:
            mat U, V;
            vec s;
            svd(U,s,V,X);
            
            mat D(ITMsize, no_basis_funcs, fill::zeros);
            mat Sig(ITMsize, ITMsize, fill::zeros);

            for(int i=0; i<s.n_rows; ++i){
                D(i,i) = s(i);
                Sig(i,i) = 1.0/(s(i)*s(i) + lambda*lambda );
            }
            
            beta.col(t) = V * D.t() * Sig * U.t() * vec(ITMY) ;
        } // Backward loop
    } // Scope
    print_DEBUG("Presim backward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    begin_time = clock();
    
    // Main-simulation
    vector<T> F;

    T dt( 1.0/double(nSteps_y) );
    T sqDt = sqrt(dt);

    // Main-simulation
    T end_result(0.0);
    
    vector<double> swap_rates(nPaths);
    mrg32k3a myRNG(seed2, seed1);
    myRNG.init(dim_n);
    vector<double> gauss(dim_n);
    //myRNG.init(M);
    //vector<double> gauss(M);
    
    number::tape->mark();
    for(int i = 0; i<nPaths; ++i){
        number::tape->rewindToMark();
        double eta = 1.0;
        T swap_path(0.0);
        for(int t = 1; t<exTimes.size(); ++t ){
            //print("exTime is ", exTimes[t], " t is ", t);
            int int_ExTime = (int)(exTimes[t]);
            int nSteps = nSteps_y * double(exTimes[t] - exTimes[t-1]);
            T swap_val(0.0);
            vector<T> lnF = log(initF);
            for(int n = 0; n<nSteps; ++n) {
                myRNG.nextG(gauss);
                //vector<T> gauss_n = MatVecProd( B_num, gauss);
                vector<T> gauss_corr = MatVecProd( B_num, gauss);
                //vector<T> gauss_cov = MatVecProd(lower, gauss);
                for(int k = int_ExTime; k < int_Tb; ++k)
                {
                    T sum(0.0);
                    for(int j = int_ExTime; j <= k; ++j ) {
                        T Fj(exp(lnF[j]));
                        sum += corr_vol[k][j] * Fj / (1.0 + tau * Fj);
                    }
                    //lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*gauss_n[k-int_Ta];
                    lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*vol[k][k]*gauss_corr[k-int_Ta];
                    //lnF[k] += vol[k][k]*tau*sum*dt - vol[k][k]*vol[k][k]/2.0*dt + sqDt*gauss_cov[k-int_Ta];
                    //print(lnF[k].value());
                } // rates
            } // steps
            // Now have F_k(T_alpha) for k=10,..,M
            F = exp(lnF);

            T floating_swap_rate;
            floating_swap_rate = SR_from_F(F, yearly_payments, int_ExTime, int_Tb );

            swap_val = disc * notional *
                C_ab( F, yearly_payments, int_ExTime, int_ExTime, int_Tb) *
                max(floating_swap_rate - r_fix, 0.0);

            swap_path += eta * swap_val;
            
            double EY2 = t < exTimes.size() - 1 ?
                // continuation value - at last ex time the eta value computed below will never be used (so 0.0 could be anything), just avoid beta of last t
                as_scalar(vec( { 1.0, double(floating_swap_rate), double(F[ int_ExTime ]) }).t() * beta.col(t)) : 0.0;
            eta *= EY2 > swap_val ? 1.0 : 0.0; // 1 if continue, 0 if exercise
            //print("Swap value ", swap_val.value());
        }
        //print("Swap path is ", swap_path.value());
        //Chol2Adj(cov_s, B_num, lower);
        
        /*for(int i = 0; i<M; i++){ // Rows
            for(int j = 0; j<dim_n; j++){ // cols
                lower[i][j].adjoint() = 0.0;
            }
        }
        vector<vector<T>> lower_test = lower;
        
        for(int i = 0; i<M; i++){ // Rows
            for(int j = 0; j<dim_n; j++){ // cols
                lower_test[i][j].adjoint() = lower[i][j].adjoint()/double(i+1);
            }
        }*/
        
        
        swap_path.propagateToMark();
        //Chol2Adj(cov_s, B_num, lower);

        end_result += swap_path;
        //print("swap_path is ", double(swap_path));
    }
    print_DEBUG("Main forward took ", float(clock() - begin_time) / CLOCKS_PER_SEC, " to compute");
    print("Simulation done ");
    
    print("lower");
    print_adj(lower);
    print("cov_s");
    print_adj(cov_s);
    print("B_num");
    print_adj(B_num);
    
    //Chol2Adj(cov_s, B_num, lower);
    print("Chol2Adj");
    
    print("lower");
    print_adj(lower);
    print("cov_s");
    print_adj(cov_s);
    
    vector<vector<T>> cov_test(M, vector<T>(M, T(0.0)));
    
    for(int i = 0; i<M; i++){ // Rows
        for(int j = 0; j<=i; j++){ // cols
            cov_test[i][j] = cov_s[i][j].adjoint();
        }
    }
    
    vector<vector<T>> G = MatMatProd(cov_test, transp(cov_test));
    print("G");
    print_val(G);
    
    mat Pert_m(Pert.size(), Pert[0].size()), G_m(G.size(), G[0].size()) ;
    
    for(int i = 0; i<Pert.size(); ++i){
        for(int j = 0; j<Pert[0].size(); ++j){
            Pert_m(i,j) = double(Pert[i][j]);
            G_m(i,j) = double(G[i][j]);
        }
    }
    
    mat U1;
    mat V1;
    
    mat U_L, V_L;
    vec s_L;
    svd(U_L,s_L,V_L,Cov);
    
    U1 = U_L.head_cols(dim_n);
    
    mat A_bar = 0.5 * (Pert_m.t()*G_m*Pert_m - ( eye(size(G_m)) - U1*U1.t() ) *
                       Pert_m.t()*G_m*Pert_m * ( eye(size(G_m)) - U1*U1.t() ) ) ;
    
    print(A_bar);
    mat PAPt = Pert_m * A_bar * Pert_m.t();
    print(PAPt );
    
    for(int i = 0; i<M; i++){ // Rows
        for(int j = 0; j<M; j++){ // cols
            cov_s[i][j].adjoint() += PAPt(i,j);
        }
    }
    

    
    number::propagateMarkToStart();
    print("PropagatedMarkToStart");
    
    print("lower");
    print_adj(lower);
    print("cov_s");
    print_adj(cov_s);

    
    //print("Corr pre prop");
    //print_adj(corr);
    for(int i = int_Ta; i<int_Tb; ++i){
        for(int j = int_Ta; j<int_Tb; ++j){
            int i2 = i-int_Ta;
            int j2 = j-int_Ta;
            
            
            //print_adj(cov_s[i2][j2]);
            //print( prev(cov_s[i2][j2].getNode())->adjoint() );
            // Propagate cov_s node
            /*cov_s[i2][j2].getNode()->propagateOne();
            
            //print_adj(cov_s[i2][j2]);
            //print( prev(cov_s[i2][j2].getNode())->adjoint() );
            // Propagate previous node
            auto prevNode = cov_s[i2][j2].getNode();
            for(int k = 0; k<10; k++){
                prevNode = prev(prevNode);
                prevNode->propagateOne();
            }*/
        }
    }
    //print("Corr post prop");
    //print_adj(corr);
    
    return end_result/double(nPaths) ;
}

#endif /* BermudaAAD_h */
