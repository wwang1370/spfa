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
