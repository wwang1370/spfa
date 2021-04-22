/* Constructor and Public Methods for Class Item
 *
 * Author: Yang Liu
 *
 * Last modified: 04/17/2021 */

#include "item.h"

/* constructor */ 

Item::Item(
  const vec &dat_,  // data vector (vec, dim = n_shortpar)
  double na_,  // missing code (double)
  const uword type,  // continuous = 0, discrete = 1 (int)
  const vec &shortpar_,  // starting parameter values (vec, dim = n_shortpar)
  const uvec &pos_,  // positivity constraint (uvec, dim = n_shortpar) 
  const uword dim_,  // dimension indicator (int)
  Basis &basis_x_,  // basis for x (Basis)
  const mat &trans_x_,  // trans matrix from shortpar to par, x (mat)
  const mat &pen_x_,  // penalty matrix, x (mat)
  const Quad &quad_x_,  // quadrature for x (Quad)
  mat &estep_wt_,  // E-step posterior weights (mat&, dim = n_obsn x n_quad)
  uword maxit_start  // maximum number of starting value iterations (int)
  ) : 
  dat(dat_), na(na_), shortpar(shortpar_), pos(pos_), dim(dim_),
  basis_x(basis_x_), trans_x(trans_x_), pen_x(pen_x_), 
  quad_x(quad_x_), estep_wt(estep_wt_)
{
  // basis, trans matrices, and quadrature for y
  switch(type)
  {
    case 0:  // continuous
      basis_y = new Bspline(basis_x.n_basis, 4, 0.0, 1.0);  // B-spline basis
      quad_y = new GaussLegendre(quad_x.n_quad1, 1, 0.0, 1.0);  // GL quadrature
      trans_y = trans_x; // transformation matrix
      pen_y = pen_x;  // penalty matrix
      break;
    case 1:  // dichotomous
      basis_y = new Iden(2);  // B-spline basis
      quad_y = new Rect(2);  // (unnormalized) rectangular quadrature
      trans_y = zeros(2, 1); // transformation matrix
      trans_y(1, 0) = 1.0;
      pen_y = zeros(1, 1);  // penalty matrix
      break;
    default:
      throw runtime_error("Item type not supported.");
      break;
  }

  // dimensions
  n_par_x = trans_x.n_rows;
  n_par_y = trans_y.n_rows;
  n_shortpar_x = trans_x.n_cols;
  n_shortpar_y = trans_y.n_cols;
  n_shortpar = n_shortpar_y + n_shortpar_y * n_shortpar_x;
  n_par = n_par_y + n_par_y * n_par_x;
   
  // gradient, search direction, and minus Hessian
  grad.set_size(n_shortpar);
  dir.zeros(n_shortpar);
  hess.set_size(n_shortpar, n_shortpar);

  // log normalizing constant
  log_norm_const.set_size(quad_x.n_quad);

  // positivity constraint
  if ( all(pos == 0) )  // no positivity constraint
  {
    search_dir_ptr = &Item::search_dir0;
  }
  else
  {
    search_dir_ptr = &Item::search_dir1;
    // active set
    activ.set_size(n_shortpar); activ.fill(2);
    activ.elem( find(pos != 0) ) = shortpar.elem( find(pos != 0) ) < 0.0;
    shortpar.elem( find(activ == 1) ).fill(0.0);
  }
}

/* copy constructor */ 

Item::Item(
  const Item &item  // reference to an Item object
  ) : 
  dat(item.dat), na(item.na), shortpar(item.shortpar), pos(item.pos), dim(item.dim),
  basis_x(item.basis_x), trans_x(item.trans_x), pen_x(item.pen_x), 
  quad_x(item.quad_x), estep_wt(item.estep_wt)
{
  basis_y = item.basis_y->clone();  // deep copy basis_y
  trans_y = item.trans_y;
  pen_y = item.pen_y;
  quad_y = item.quad_y->clone();  // deep copy quad_y
  n_par = item.n_par; 
  n_par_y = item.n_par_y;
  n_par_x = item.n_par_x;
  n_shortpar = item.n_shortpar;
  n_shortpar_y = item.n_shortpar_y;
  n_shortpar_x = item.n_shortpar_x;
  par = item.par;
  f = item.f;
  grad = item.grad;
  hess = item.hess;
  dir = item.dir;
  pen_val = item.pen_val;
  log_norm_const = item.log_norm_const;
  activ = item.activ;
  search_dir_ptr = item.search_dir_ptr;
}

/* destructor */

Item::~Item()
{
  delete basis_y;
  delete quad_y;
}

/* extend_par: shortpar -> par
 *
 * updates: par (double, dim = n_par) */

void Item::extend_par()
{
  par.set_size(n_par);
  mat spmat(shortpar.memptr(), n_shortpar_y, n_shortpar_x + 1, false),
    pmat(par.memptr(), n_par_y, n_par_x + 1, false);
  pmat.col(0) = trans_y * spmat.col(0);
  pmat.tail_cols(n_par_x) = trans_y * 
    spmat.tail_cols(n_shortpar_x) * trans_x.t();
}

/* reduce_par: par -> shortpar
 *
 * updates: shortpar (double, dim = n_shortpar) */

void Item::reduce_par()
{
  shortpar.set_size(n_shortpar);
  mat spmat(shortpar.memptr(), n_shortpar_y, n_shortpar_x + 1, false),
    pmat(par.memptr(), n_par_y, n_par_x + 1, false);
  mat pinv_y = pinv(trans_y), pinv_x = pinv(trans_x);   // compute pseudoinverses
  spmat.col(0) = pinv_y * pmat.col(0);
  spmat.tail_cols(n_shortpar_x) = pinv_y * 
    pmat.tail_cols(n_par_x) * pinv_x.t();
}
