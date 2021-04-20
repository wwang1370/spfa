/* Wrapper functions
 *
 * Author: Yang Liu
 *
 * Last modified: 04/13/2021 */

#include "test.h"

/* spfa_main: fitting semi-parametric factor model
 *
 * returns: stuffs generated by test.output (list) */

// [[Rcpp::export]]
Rcpp::List spfa_main(
  const arma::mat &dat,   // data matrix (mat, dim = n_obsn x n_item)
  double na,  // missing code (double) 
  const arma::uvec &item_type,  // item type (uvec, dim = n_item)
  const Rcpp::List &shortpar,  // starting values (mat, dim = max(n_shortpar) x n_item)
  const Rcpp::List &pos,  // positivity constraints (mat, dim = max(n_shortpar) x n_item)  
  arma::uword n_basis,  // number of basis functions (int)
  double lmbd,  // penalty weights for y (double)
  arma::uword n_quad,  // number of quadrature points (int) */
  const arma::uvec &dim,  // dimension indicator (uvec, dim = n_item)
  arma::uword maxit_em,  // max number of EM iterations (int)
  arma::uword maxit_mstep,  // max number of M-step iterations (int)
  arma::uword maxit_start,  // max number of starting value iterations (int)
  double tol_em,  // tolerance for likelihood change, EM (double) 
  double tol_mstep,  // tolerance for gradient norm, M-step (double) 
  double tol_start,  // tolerance for gradient norm, starting values (double) 
  int n_thrd  // number of threads (int)
  )
{
  // initialization
  Test test(dat, na, item_type, shortpar, pos, 
    n_basis, lmbd, 
    n_quad, recode_unique(dim),
    maxit_em, maxit_mstep, maxit_start, 
    tol_em, tol_mstep, tol_start, 
    n_thrd);
  // EM 
  test.em();
  // output
  Rcpp::List ret = test.output();
  return ret;
}

///* spfa_score: scoring based on the semi-parametric model
// *
// * returns: latent variable scores (double, dim = n_obsn) */
//
//// [[Rcpp::export]]
//arma::mat spfa_score(
//  const arma::mat& shortpar,  // starting values (double&, dim = n_shortpar x n_item)
//  const arma::mat& dat,   // data matrix (int&, dim = n_obsn x n_item)
//  arma::uword n_basis,  // number of basis functions (int)
//  arma::uword n_quad,  // number of quadrature points (int) */
//  arma::uword mode  // mode (uint): 0 = uniform, 1 = normal
//  )
//{
//  // initialization
//  uvec type(shortpar.n_cols); type.fill(0);
//  Test test(shortpar, dat, type, n_basis, 0.0, n_quad, 
//    0, 0, 0, 0.0, 0.0, 0.0);
//  test.start_val();  // starting values (only compute log_norm_const)
//  test.estep(); // run E-step to get weights
//  arma::mat x = test.score(mode);  // compute scores
//  return x;
//}
//
///* marg_lik: compute marginal likelihood for each observation
// *
// * returns: marginal likelihood (vec, dim = n_obsn) */
//
//// [[Rcpp::export]]
//arma::vec marg_lik(
//  const arma::mat& shortpar,  // starting values (double&, dim = n_shortpar x n_item)
//  const arma::mat& dat,   // data matrix (int&, dim = n_obsn x n_item)
//  arma::uword n_basis,  // number of basis functions (int)
//  arma::uword n_quad  // number of quadrature points (int) */
//  )
//{
//  // test initialization
//  arma::uword n_item = shortpar.n_cols;
//  uvec type(n_item); type.fill(0);
//  Test test(shortpar, dat, type, n_basis, 0.0, n_quad, 
//    0, 0, 0, 0.0, 0.0, 0.0);
//  test.start_val();  // starting values (only compute log_norm_const)
//  test.estep(); // run E-step to get weights
//  arma::uvec it = arma::regspace<uvec>(0, n_item - 1);
//  arma::vec f = test.marg_lik(dat, it);
//  return f;
//}
//
/////* cond_dns: evaluate conditional density
//// *
//// * returns: deviance (double) */
////
////// [[Rcpp::export]]
////arma::mat cond_dns(
////  arma::vec shortpar,  // parameter values (vec, dim = n_shortpar)
////  arma::vec y,  // y values (double, dim = y.n_elem)
////  arma::vec x,  // x values (double, dim = x.n_elem)
////  arma::uword n_basis,  // number of basis functions (int)
////  arma::uword n_quad  // number of quadrature points (int)
////  )
////{
////  Bspline bspl(n_basis, 4, 0.5, 0.0);
////  Quad quad(n_quad, 0.0, 1.0);
////  arma::mat estep_wt;
////  Item item(shortpar, y.memptr(), 0, 0, bspl, quad, estep_wt, 0);
////  arma::mat f = trunc_exp( item.cond_log_dns(y, x) );
////  return f;
////}
//
/* cubic_bspl: cubic B-spline basis matrix
 * 
 * returns: basis matrix (mat, dim = n x n_basis) */

// [[Rcpp::export]]
arma::mat cubic_bspl(
  arma::vec x,   // data vector (int, dim = n)
  arma::uword n_basis  // number of basis functions (int)
  )
{
  Bspline basis(n_basis, 4, 0.0, 1.0);
  arma::mat ret(x.n_elem, n_basis);
  for (arma::uword i = 0; i < x.n_elem; ++i)
    ret.row(i) = basis.eval( x(i) );
  return ret;
}
//
///* cubic_bspl0: cubic B-spline basis matrix, with side condition
// * 
// * returns: basis matrix (mat, dim = n x (n_basis - 1) ) */
//
//// [[Rcpp::export]]
//arma::mat cubic_bspl0(
//  arma::vec x,   // data vector (int, dim = n)
//  arma::uword n_basis,  // number of basis functions (int)
//  double x0 // x value where side condition applies (double)
//  )
//{
//  Bspline basis(n_basis, 4, x0, 1);
//  arma::mat ret(x.n_elem, n_basis - 1);
//  for (arma::uword i = 0; i < x.n_elem; ++i)
//    ret.row(i) = basis.eval0( x(i) );
//  return ret;
//}

/* gl_quad: set up the rescaled gauss-Legendre quadrature
 * 
 * returns: nodes and weigths (list) */

// [[Rcpp::export]]
Rcpp::List gl_quad(
  arma::uword n_quad,  // number of quadarture points (int)
  arma::uword dim,  // dimension (int)
  double lower,  // lower limit (double)
  double upper  // upper limit (double)
  )
{
  GaussLegendre gl(n_quad, dim, lower, upper);
  List ret = List::create(
    Named("x") = gl.node,
    Named("w") = gl.weight 
    );
  return ret;
}

/* reduce_par: par -> shortpar
 * extend_par: shortpar -> par
 * 
 * returns: transformed parameters (double) */

// [[Rcpp::export]]
arma::vec reduce_par(
  arma::uword n_basis,  // number of basis functions (int)
  arma::vec par  //  parameters (double, length = n_par)
  )
{
  // prepare quantities
  arma::vec dat; arma::uvec pos; arma::mat estep_wt;
  Bspline bspl(n_basis, 4, 0.0, 1.0);
  arma::rowvec b0 = bspl.eval(0.5);
  mat tr = null(b0).t();
  tr = solve(tr * diff_mat(n_basis, 1).t(), tr);
  inplace_trans(tr);
  arma::mat pen(tr.n_cols, tr.n_cols);
  GaussLegendre quad(2, 1, 0.0, 1.0);
  // initialize item
  Item item(dat, -9.0, 0, {0}, pos, 0,
    bspl, tr, pen, quad, estep_wt, 0);
  // call public methods
  item.set_par(par);
  item.reduce_par();
  return item.get_shortpar();
}

// [[Rcpp::export]]
arma::vec extend_par(
  arma::uword n_basis,  // number of basis functions (int)
  arma::vec shortpar  //  short parameters (double, length = n_shortpar)
  )
{
  // prepare quantities
  arma::vec dat; arma::uvec pos; arma::mat estep_wt;
  Bspline bspl(n_basis, 4, 0.0, 1.0);
  arma::rowvec b0 = bspl.eval(0.5);
  mat tr = null(b0).t();
  tr = solve(tr * diff_mat(n_basis, 1).t(), tr);
  inplace_trans(tr);
  arma::mat pen(tr.n_cols, tr.n_cols);
  GaussLegendre quad(2, 1, 0.0, 1.0);
  // initialize item
  Item item(dat, -9.0, 0, shortpar, pos, 0,
    bspl, tr, pen, quad, estep_wt, 0);
  // call public methods
  item.extend_par();
  return item.get_par();
}
