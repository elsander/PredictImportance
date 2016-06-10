#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Simulate a discrete-time Lotka-Volterra time series in C
//'
//' Simulate discrete-time Lotka-Volterra population dynamics
//' with perturbed growth rates
//'
//' @param r_ a matrix with S (number of species) rows and simtime
//' cols containing the growth rates at each time step (potentially
//' perturbed at each time step) 
//' @param alphas_ an S by S matrix of per capita interaction strengths
//' @param n0s_ a vector of length S containing species abundances
//' @param Is_ a vector of length S containing immigration rates
//' @param deltat the size of the time step
//' @param simtime the total number of time steps to simulate
//'
//' @return an S by simtime matrix of species abundances at each timestep
//' @export
// [[Rcpp::export]]
NumericMatrix discreteLV_C(NumericMatrix r_,
                           NumericMatrix alphas_,
                           NumericVector n0s_,
                           NumericVector Is_,
                           double deltat,
                           int simtime){
    // copy the data into armadillo structures
    arma::mat r = Rcpp::as<arma::mat> (r_);
    arma::mat alphas = Rcpp::as<arma::mat>(alphas_);
    arma::colvec n0s = Rcpp::as<arma::colvec>(n0s_);
    arma::colvec Is = Rcpp::as<arma::colvec>(Is_);

    int S = n0s.size();
    NumericMatrix ns(S, simtime);
    int i,j,k;

    // used for matrix calculations
    arma::colvec deltatvec(S);
    arma::colvec tmp(S);
    arma::colvec Imm_n0s(S);

    for(i=0; i<S; i++){
        deltatvec(i) = deltat;
    }

    for(i=0; i<(simtime-1); i++){
        // copy current state of n into output matrix ns
        for(k=0; k<S; k++){
            ns(k,i) = n0s(k);
        }
        // check for user interrupts
        Rcpp::checkUserInterrupt();
        // The number of intermediate steps is determined by deltat
        for(j=0; j<(int)(1/deltat); j++){
            
            // calculate next abundance for each species
            // The % operator is for elementwise multiplication
            tmp = deltatvec % (r.col(i) + (alphas * n0s) + (Is / n0s));
            for(k=0; k<S; k++){
                tmp(k) = exp(tmp(k));
            }
            n0s = n0s % tmp;
        }
    }

    // final timepoint
    for(i=0; i<S; i++){
        ns(i,simtime-1) = n0s(i);
    }
    
    return ns;
}
