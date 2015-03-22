#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
MatrixXd sherman(const Map<MatrixXd> Ap, const Map<MatrixXd> u, const Map<MatrixXd> v) {
  // R formula for reference:
  //Ap - (Ap %*% u %*% t(v) %*% Ap)/drop(1 + t(v) %*% Ap %*% u)
  MatrixXd a(Ap.rows(), Ap.cols());
  double b;
  MatrixXd vt(1, Ap.rows());
  vt = v.transpose();

  return Ap - (Ap * u * vt * Ap)/(1 + (vt * Ap * u)(0, 0));
}

// [[Rcpp::export]]
MatrixXd doublesherman(const Map<MatrixXd> Ap, const Map<MatrixXd> u, const Map<MatrixXd> v) {
  // R formula for reference:
  //Sherman is solve(A + u %*% t(v) ) =  Ap - (Ap %*% u %*% t(v) %*% Ap)/drop(1 + t(v) %*% Ap %*% u)
  MatrixXd a(Ap.rows(), Ap.cols());
  double b;
  MatrixXd vt(1, Ap.rows());
  MatrixXd shermanterm(Ap.rows(), Ap.cols()) ;

  vt = v.transpose();

  shermanterm = (Ap * u * vt * Ap)/(1 + (vt * Ap * u)(0, 0)) ;
  return Ap - shermanterm - shermanterm.transpose() ;
}
// There's more savings to be had here, when v has only a single non-zero entry 
