// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// [[Rcpp::export]]
arma::mat XbetaZb_cpp(const arma::mat &X, const arma::vec &beta, 
            const arma::mat &Z, const arma::vec &b) {
  return (X * beta + Z * b);
}

// [[Rcpp::export]]
arma::mat Xbeta_cpp(const arma::mat &X, const arma::vec &beta) {
  return (X * beta);
}

// [[Rcpp::export]]
double dnorm_cpp(const arma::vec &y, const arma::mat &X, const arma::vec &beta, 
                 const arma::mat &Z, const arma::vec &b, double sigma){
    
    double variance = sigma * sigma;
    double expPart = sum(square(y - X * beta - Z * b));
    
    return (-0.5 * (y.n_elem * log(2 * M_PI) + y.n_elem * log(variance) + expPart/variance));
}

// [[Rcpp::export]]
arma::vec dmvnorm_cpp(const arma::mat& x, const arma::rowvec& mean, 
                   const arma::mat& sigma, bool logd = false) {
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double log2pi = std::log(2.0 * M_PI);
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    for (int i = 0; i < n; ++i) {
        arma::vec z = rooti * arma::trans(x.row(i) - mean) ;
        out(i) = constants - 0.5 * arma::sum(z%z) + rootisum;
    }
    if (logd == false) {
        out = exp(out);
    }
    return(out);
}

// [[Rcpp::export]]
double logSurvival_cpp(const arma::mat &Xs, const arma::vec &betas_new, const arma::mat &Zs, 
                       const arma::vec &b, const arma::vec &alphas_new, const arma::vec &W_time_independent, 
                       const arma::vec &gammas_new, const arma::mat &W2s, const arma::vec &Bs_gammas_new, 
                       const arma::vec &wk,  double P) {
  
  arma::mat Ys, Ys_extra;
  arma::vec longitudinal_part;
  arma::vec Vi, logbaseline_hazard;
  
  Ys = Xs * betas_new + Zs * b;
  longitudinal_part = Ys * alphas_new;
  
  logbaseline_hazard = W2s * Bs_gammas_new;
  Vi = exp(logbaseline_hazard + longitudinal_part);
  
  double exp_eta_tw = exp(sum(W_time_independent % gammas_new));
  
  return exp_eta_tw * P * sum(wk % Vi);
}
