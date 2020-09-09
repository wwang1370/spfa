/* Constructors: Class _Item_
 *
 * Author: Yang Liu
 *
 * Last modified: 08/27/2020 */

#include "item.h"

/* constructor */ 

Item::Item(
  const vec& shortpar_,  // starting parameter values (vec, dim = n_shortpar)
  double* dat_ptr,  // data vector (vec&)
  uword n_obsn,  // number of observations (int)
  const uword type_,  // item type (int, dim = 1) 
  Bspline& bspl_,  // B-spline basis (Bspline&) 
  Quad& quad_,  // quadrature (Bspline&)
  mat& estep_wt_,  // E-step posterior weights (mat&, dim = n_obsn x n_quad)
  uword maxit_start  // maximum number of starting value iterations (int)
  ) : 
  type(type_), 
  bspl(bspl_), 
  quad(quad_),
  estep_wt(estep_wt_)
{
  // dimensions
  n_par = bspl.n_basis * (1 + bspl.n_basis);
  n_shortpar = (bspl.n_basis - 1) * bspl.n_basis;
  // data
  // workaround to the segfault error (YL 08/31/2020)
  dat = vec(dat_ptr, n_obsn, false);
  // parameters, gradient, search direction, and minus Hessian
  shortpar = shortpar_;
  grad.set_size(n_shortpar);
  dir.zeros(n_shortpar);
  hess.set_size(n_shortpar, n_shortpar);
  // log normalizing constant
  log_norm_const.set_size(quad.n_quad);
  // function pointers
  switch (type)
  {
    case 0:  // type 0: no constraint
      search_dir_ptr = &Item::search_dir0;
      break;
    case 1:  // type 1: likelihood-ratio ordering
      search_dir_ptr = &Item::search_dir1;
      // active set
      activ.set_size(n_shortpar); activ.fill(2);
      activ.subvec(bspl.n_basis - 1, n_shortpar - 1) =
        shortpar.subvec(bspl.n_basis - 1, n_shortpar - 1) < 0.0;
      shortpar.elem( find(activ == 1) ).fill(0.0);
      break;
    default: 
      Rcout << "Error: Unsupported item type " << type << endl;
      exit(1);
      break;
  }
}
