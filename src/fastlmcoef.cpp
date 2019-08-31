#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
//[[Rcpp::export]]




arma::vec fastLMcoef(arma::mat X, arma::colvec y){
  arma::vec coef = arma::solve(X, y);
  return(coef);
  //only returns the coefficients, taken from Dirk Eddelbuettel
}
