/* Utility functions
 *
 * Author: Yang Liu
 *
 * Last modified: 09/30/2020 */

#include "util.h"

/* choose: compute combinatory number 
 *
 * return: combinatory number (int) */

uword choose(
  uword n,   // total size (int)
  uword k  // subset size (int)
  )
{
  double n_ = (double) n, k_ = (double) k;
  double y = lgamma(n_ + 1.0) - lgamma(k_ + 1.0) - lgamma(n_ - k_ + 1.0);
  uword ret = (uword) (std::exp(y) + 0.5);
  return ret;
}

/* grid_loc: convert index to grid locations
 *
 * return: grid location (uvec, dim = n_dim) */

uvec grid_loc(
    uword indx,  // index to be converted to grid location (int)
    uword n_dim,   // number of dimensions (int)
    uword n_pts   // number of points per dimension (int)
  )
{
  uvec loc(n_dim);
  for (uword d = 0; d < n_dim; ++d)
  {
    loc(d) = indx % n_pts;
    indx = ( indx - loc(d) ) / n_pts;
  }
  return loc;
}

/* pow_uword: raise unsigned integer to an unsigned integer power 
 *
 * return: result (uword) */

uword pow_uword(
  uword base,  // base number (uword)
  uword exp  // exponent (uword)
  )
{
  uword ret = base;
  for (uword d = 1; d < exp; ++d)
    ret *= base;
  return ret;
}

/* expand_grid: expand vector x to an n_dim-dimensional grid
 *
 * return: grid location (matrix, dim = n_pts^n_dim x n_dim) */

mat expand_grid(
  vec x,  // vector to be expanded (vec, dim = n_pts)
  uword n_dim   // number of dimensions (int)
  )
{
  uword n_pts = x.n_elem, n_grid = pow_uword(n_pts, n_dim);
  mat ret(n_grid, n_dim);
  for (uword i = 0; i < n_grid; ++i)
  {
    uvec grid_indx = grid_loc(i, n_dim, n_pts);
    for (uword d = 0; d < n_dim; ++d)
      ret(i, d) = x( grid_indx(d) );
  }
  return ret;
}

/* diff_mat: evaluate the difference matrix of a given order
 *
 * returns: differencing matrix (double, dim = (n - order) x n) */

mat diff_mat(
  uword n,  // dimension (int)
  uword order  // order of difference (int)
  )
{
  mat diff(n - order, n); diff.zeros();
  if (order == 1)  // exit: order = 1
  {
    diff.diag(0).fill(1.0);
    diff.diag(1).fill(-1.0);
  }
  else  // otherwise, do recursion
  {
    // apply 1st-order diff matrix to columns of diff_mat(n, order - 1)
    mat diff_ = diff_mat(n, order - 1);
    for (uword j = 0; j < n; ++j)
    {
      uword i = j < order ? 0 : j - order;
      while (i < j + 1 && i < diff.n_rows)
      {
        diff(i, j) = diff_(i, j) - diff_(i + 1, j);
        ++i;
      }
    }
  }
  return diff;
}

/* is_equal: float number comparison
 *
 * return: equal or not (bool) */

bool is_equal(
  double a,  // first scalar (double)
  double b  // second scalar (double)
  )
{
  if (std::abs(a - b) < DBL_EPS) return true;
  else return false;
}

/* arma2r: convert arma::vec to Rcpp vector
 *
 * return: R vector (NumericVector) */

NumericVector arma2r(
  vec x  // armadillo vector (vec)
  )
{
  NumericVector y = Rcpp::wrap(x);
  y.attr("dim") = R_NilValue;
  return y;
}

/* recode_unique: recode sequence x to 0, ..., unique(x).n_elem - 1
 *
 * return: recoded vector (uvec) */

uvec recode_unique(
  uvec x   // unsigned integer sequence to be recoded (uvec)
  )
{
  uvec u = unique(x), y( size(x) );
  for (uword i = 0; i < u.n_elem; ++i)
    y.elem( find( x == u(i) ) ).fill(i);
  return y;
}
