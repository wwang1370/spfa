/* Wrapper functions
 *
 * Author: Yang Liu
 *
 * Last modified: 05/19/2020 */

#include "test.h"

/* spfa_main: fitting the semi-parametric factor model */

// [[Rcpp::export]]
Rcpp::List spfa_main(
  const arma::mat& shortpar,  // starting values (double&, dim = n_shortpar x m)
  arma::mat& dat,   // data matrix (int&, dim = n x m)
  const arma::uvec& type,  //item type (int, length = m)  
  arma::uword n_basis,  // number of basis functions (int)
  double lmbd,  // penalty weights for y (double)
  arma::uword n_quad,  // number of quadrature points (int) */
  arma::uword maxit_em,  // max number of EM iterations (int)
  arma::uword maxit_mstep,  // max number of M-step iterations (int)
  arma::uword maxit_start,  // max number of starting value iterations (int)
  double tol_em,  // tolerance for likelihood change, EM (double) 
  double tol_mstep,  // tolerance for gradient norm, M-step (double) 
  double tol_start  // tolerance for gradient norm, starting values (double) 
  )
{
  // initialization
  Test test(shortpar, dat, type, n_basis, lmbd, n_quad, 
    maxit_em, maxit_mstep, maxit_start, tol_em, tol_mstep, tol_start);
  // starting values
  test.start_val();
  // EM 
  test.em();
  // output
  Rcpp::List ret = test.output();
  return ret;
}

/* cubic_bspl: cubic B-spline basis matrix
 * 
 * returns: basis matrix (mat, dim = n x n_basis) */

// [[Rcpp::export]]
arma::mat cubic_bspl(
  arma::vec x,   // data vector (int, dim = n)
  arma::uword n_basis  // number of basis functions (int)
  )
{
  Bspline basis(n_basis, 4);
  arma::mat ret = basis.eval(x);
  return ret;
}

/* cubic_bspl0: cubic B-spline basis matrix, with side condition
 * 
 * returns: basis matrix (mat, dim = n x (n_basis - 1) ) */

// [[Rcpp::export]]
arma::mat cubic_bspl0(
  arma::vec x,   // data vector (int, dim = n)
  arma::uword n_basis,  // number of basis functions (int)
  double x0 // x value where side condition applies (double)
  )
{
  Bspline basis(n_basis, 4, x0, 1);
  arma::mat ret = basis.eval0(x);
  return ret;
}

/* gl_quad: set up the rescaled gauss-Legendre quadrature
 * 
 * returns: nodes and weigths (list) */

// [[Rcpp::export]]
Rcpp::List gl_quad(
  arma::uword n_quad,  // number of quadarture points (int)
  double lower,  // lower limit (double)
  double upper  // upper limit (double)
  )
{
  Quad quad(n_quad, lower, upper);
  List ret = List::create(
    Named("x") = quad.node,
    Named("w") = quad.weight);
  return ret;
}
