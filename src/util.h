/* Headers: Utility functions
 *
 * Author: Yang Liu
 *
 * Last modified: 09/30/2020 */

#ifndef UTIL_H
#define UTIL_H

#include "RcppArmadillo.h"
using namespace std;
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

uword fact(uword n);  // compute factorial
uword choose(uword n, uword k);  // compute combinatory number
uvec grid_loc( uword indx, uword n_dim, uword n_pts); // index -> grid locations

#endif
