/* Headers: Utility functions
 *
 * Author: Yang Liu
 *
 * Last modified: 04/13/2021 */

#ifndef UTIL_H
#define UTIL_H

#include "RcppArmadillo.h"
#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
  #define omp_get_max_threads() 1
  #define omp_get_num_threads() 1
#endif

using namespace std;
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// constants
static const double DBL_EPS = 1e-64;   // tolerance for comparison between floats

uword fact(uword n);  // compute factorial
uword choose(uword n, uword k);  // compute combinatory number
uvec grid_loc(uword indx, uword n_dim, uword n_pts); // index -> grid locations
uword pow_uword(uword base, uword exp); // function to compute integer power
mat expand_grid(vec x, uword n_dim); // expand x to a grid of dimension n_dim
mat diff_mat(uword n, uword order);  // differencing matrix
bool is_equal(double a, double b);  // float number comparison
uvec recode_unique(uvec x);  // recode sequence x to 0, ..., unique(x).n_elem - 1
NumericVector arma2r(vec x);  // arma::vec to Rcpp::NumericVector

#endif
