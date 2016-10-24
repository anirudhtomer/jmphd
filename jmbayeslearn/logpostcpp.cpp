// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double temp_cpp(const arma::mat &Xs, const arma::vec &betas_new, const arma::mat &Zs, 
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
