/* Wrapper functions
 *
 * Author: Yang Liu
 *
 * Last modified: 09/30/2020 */

#include "test.h"

/* spfa_main: fitting the semi-parametric factor model 
 *
 * returns: stuffs generated by test.output (list) */

// [[Rcpp::export]]
Rcpp::List spfa_main(
  const arma::mat& shortpar,  // starting values (double&, dim = n_shortpar x n_item)
  arma::mat& dat,   // data matrix (int&, dim = n_obsn x n_item)
  const arma::uvec& type,  //item type (int, length = n_item)  
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

/* spfa_score: scoring based on the semi-parametric model
 *
 * returns: latent variable scores (double, dim = n_obsn) */

// [[Rcpp::export]]
arma::vec spfa_score(
  const arma::mat& shortpar,  // starting values (double&, dim = n_shortpar x n_item)
  arma::mat& dat,   // data matrix (int&, dim = n_obsn x n_item)
  arma::uword n_basis,  // number of basis functions (int)
  arma::uword n_quad  // number of quadrature points (int) */
  )
{
  // initialization
  uvec type(shortpar.n_cols); type.fill(0);
  Test test(shortpar, dat, type, n_basis, 0.0, n_quad, 
    0, 0, 0, 0.0, 0.0, 0.0);
  test.start_val();  // starting values (only compute log_norm_const)
  test.estep(); // run E-step to get weights
  arma::vec x = test.score();  // compute scores
  return x;
}

/* xv_risk: compute cross-validation risk
 *
 * returns: risk function for each combination (vec) */

// [[Rcpp::export]]
arma::vec xv_risk(
  const arma::mat& shortpar,  // starting values (double&, dim = n_shortpar x n_item)
  arma::mat& dat,   // data matrix (int&, dim = n_obsn x n_item)
  arma::uword n_basis,  // number of basis functions (int)
  arma::uword n_quad,  // number of quadrature points (int) */
  arma::uword order  // order of margins (int)
  )
{
  // test initialization
  arma::uword n_item = shortpar.n_cols;
  uvec type(n_item); type.fill(0);
  Test test(shortpar, dat, type, n_basis, 0.0, n_quad, 
    0, 0, 0, 0.0, 0.0, 0.0);
  test.start_val();  // starting values (only compute log_norm_const)
  test.estep(); // run E-step to get weights

  // loop over all combinations of items
  arma::vec risk( choose(n_item, order) );
  string bitmask(order, 1);  // order leading 1s
  bitmask.resize(n_item, 0);  // n_item - order trailing 0s
  arma::uword c = 0;
  do 
  {
    uvec it(order);  // item combination
    uword s = 0;
    for (arma::uword i = 0; i < n_item; ++i)
    {
      if (bitmask[i]) 
      {
        it[s] = i;
        s++;
      }
      if (s >= order) break;
    }
    // print info
    Rcout << "Compute risk for items ";
    for (arma::uword d = 0; d < order; ++d)
      Rcout << it(d) << ' ';
    Rcout << '\r';
    risk(c) = test.risk(it);  // compute risk
    c++;
  } while ( prev_permutation( bitmask.begin(), bitmask.end() ) );
  Rcout << endl;
  return risk;
}

/* cond_dns: evaluate conditional density
 *
 * returns: deviance (double) */

// [[Rcpp::export]]
arma::mat cond_dns(
  arma::vec shortpar,  // parameter values (vec, dim = n_shortpar)
  arma::vec y,  // y values (double, dim = y.n_elem)
  arma::vec x,  // x values (double, dim = x.n_elem)
  arma::uword n_basis,  // number of basis functions (int)
  arma::uword n_quad  // number of quadrature points (int)
  )
{
  Bspline bspl(n_basis, 4, 0.5, 0.0);
  Quad quad(n_quad, 0.0, 1.0);
  arma::mat estep_wt;
  Item item(shortpar, y.memptr(), 0, 0, bspl, quad, estep_wt, 0);
  arma::mat f = arma::exp( item.cond_log_dns(y, x) );
  return f;
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
  vec dat; mat estep_wt;
  Bspline bspl(n_basis, 4, 0.5, 0.0);
  Quad quad(2, 0.0, 1.0);
  arma::uword n_shortpar = n_basis * (n_basis - 1);
  Item item(par.head(n_shortpar), dat.memptr(), 0, 0, bspl, quad, estep_wt, 0);
  item.par = par;
  item.reduce_par();
  return item.shortpar;
}

// [[Rcpp::export]]
arma::vec extend_par(
  arma::uword n_basis,  // number of basis functions (int)
  arma::vec shortpar  //  short parameters (double, length = n_shortpar)
  )
{
  vec dat; mat estep_wt;
  Bspline bspl(n_basis, 4, 0.5, 0.0);
  Quad quad(2, 0.0, 1.0);
  Item item(shortpar, dat.memptr(), 0, 0, bspl, quad, estep_wt, 0);
  arma::uword n_par = n_basis * (n_basis + 1);
  item.par.set_size(n_par);
  item.extend_par();
  return item.par;
}
