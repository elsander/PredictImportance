#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix discreteLV_C(NumericMatrix r,
                           NumericMatrix alphas,
                           NumericVector ntmp,
                           double deltat,
                           int simtime){
    // THIS FUNCTION ALTERS ntmp IN PLACE!
    int S = ntmp.size();
    NumericMatrix ns(S, simtime);
    int i,z,k,m;
    double j,tmp, doublei;

    for(i=0; i<(simtime-1); i++){
        // copy current state of n into output matrix ns
        for(z=0; z<S; z++){
            ns(z,i) = ntmp(z);
        }
        // take steps by deltat
        // cast i as double for use in for loop
        doublei = (double)i;
        for(j=doublei; j<(doublei+1.0); j+=deltat){
            // calculate next abundance for each species
            for(k=0; k<S; k++){
                tmp = 0;
                for(m=0; m<S; m++){
                    tmp += alphas(k,m) * ntmp(m);
                }
                // r(k,i) because we are reusing the r value for
                // the integer timestep for each deltat step
                // in between (may add splining later)
                ntmp(k) = ntmp(k)*exp(deltat * (r(k,i) + tmp));
            }
        }
    }

    for(i=0; i<S; i++){
        ns(i,simtime-1) = ntmp(i);
    }
    
    return ns;
}

///////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericMatrix discreteLV_C_arma(NumericMatrix r_,
                                NumericMatrix alphas_,
                                NumericVector ntmp_,
                                double deltat,
                                int simtime){
    // copy the data to armadillo structures
    arma::mat r = Rcpp::as<arma::mat> (r_);
    arma::mat alphas = Rcpp::as<arma::mat>(alphas_);
    arma::colvec ntmp = Rcpp::as<arma::colvec>(ntmp_);

    int S = ntmp.size();
    NumericMatrix ns(S, simtime);
    int i,j,k;
    // double j, doublei;

    // used for matrix calculations
    arma::colvec deltatvec(S);
    arma::colvec tmp(S);

    for(i=0; i<S; i++){
        deltatvec(i) = deltat;
    }

    for(i=0; i<(simtime-1); i++){
        // copy current state of n into output matrix ns
        for(k=0; k<S; k++){
            ns(k,i) = ntmp(k);
            // if(i == 0){
            //     Rprintf("ntmp_ %f ntmp %f ns %f\n", ntmp_(k), ntmp(k), ns(k,i));
            // }
        }
        // take steps by deltat
        // cast i as double for use in for loop
        // doublei = (double)i;
        // for(j=doublei; j<(doublei+1.0); j+=deltat){
        for(j=0; j<(int)(1/deltat); j++){
            // Rprintf("j %f, doublei %f\n", j, doublei);
            
            // calculate next abundance for each species
            // The % operator is for elementwise multiplication
            tmp = deltatvec % (r.col(i) + alphas * ntmp);
            // apply a lambda to exponentiate each element
            // should have similar performance to a regular for loop
            // The ampersand is there because we MUST pass by
            // reference
            for(k=0; k<S; k++){
                tmp(k) = exp(tmp(k));
            }
            ntmp = ntmp % tmp;
        }
    }

    for(i=0; i<S; i++){
        ns(i,simtime-1) = ntmp(i);
    }
    
    return ns;
}
