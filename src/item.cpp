/* Constructor and Public Methods for Class Item
 *
 * Author: Yang Liu
 *
 * Last modified: 03/08/2022 */

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
  mat &estep_wt_  // E-step posterior weights (mat&, dim = n_obsn x n_quad)
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
    case 1:  // discrete
    {
      // find the maximum data entry
      double m = -DBL_MAX;
      for (uword i = 0; i < dat.n_elem; ++i)
      {
        if ( is_equal(dat(i), na) ) continue;
        if (dat(i) > m) m = dat(i);
      }
      uword n_cat = ( (uword) m ) + 1;  // # of categories
      basis_y = new Iden(n_cat);  // identity basis
      quad_y = new Rect(n_cat);  // (unnormalized) rectangular quadrature
      trans_y = zeros(n_cat, n_cat - 1); // transformation matrix
      trans_y.diag(-1).ones();
      pen_y = zeros(n_cat - 1, n_cat - 1);  // penalty matrix
      break;
    }
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
  if (shortpar.n_elem != n_shortpar)
    throw runtime_error("Length of shortpar is not identical to n_shortpar.");
  n_par = n_par_y + n_par_y * n_par_x;
   
  // gradient, search direction, and minus Hessian
  grad.set_size(n_shortpar);
  dir.zeros(n_shortpar);
  hess.set_size(n_shortpar, n_shortpar);
  cond1 = DBL_MAX;

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
  cond1 = item.cond1;
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

/* basis_exp: evaluate basis expansion and its derivatives
 *
 * returns: basis expansion value (double)
 *
 * updates: 
 * gr: gradient (if deriv = true, destroy first; vec, dim = n_shortpar) */

double Item::basis_exp(
  vec& gr,  // gradient vector (vec&, dim = n_shortpar)
  double y,  // y value (double)
  double x,  // x value (double)
  bool deriv  // whether derivatives evaluated (bool)
  )
{
  mat spmat(shortpar.memptr(), n_shortpar_y, n_shortpar_x + 1, false);
  rowvec design_y = basis_y->eval(y) * trans_y, 
    design_x = basis_x.eval(x) * trans_x;
  // basis expansion
  double f = dot(spmat.col(0), design_y) + 
    dot(design_y * spmat.tail_cols(n_shortpar_x), design_x);
  // gradient
  if (deriv)
  {
    gr.set_size(n_shortpar);
    mat gmat(gr.memptr(), n_shortpar_y, n_shortpar_x + 1, false);
    gmat.col(0) = trans(design_y);
    gmat.tail_cols(n_shortpar_x) = gmat.col(0) * design_x;
  }
  return f;
}

/* log_normalize: evaluate log normalizing constant and derivatives
 *
 * returns: log normalizing constant value (double)
 *
 * updates: 
 * gr: gradient (if deriv = true, destroy first; vec, dim = n_shortpar) 
 * he: Hessian (if deriv = true, destroy first; mat, 
 *     dim = n_shortpar x n_shortpar) */

double Item::log_normalize(
  vec& gr,  // gradient vector (vec&, dim = n_shortpar)
  mat& he,  // Hessian matrix (mat&, dim = n_shortpar x n_shortpar)
  double x,  // x value (double)
  bool deriv  // whether derivatives are evaluated (bool)  
  )
{
  vec udns(quad_y->n_quad);  // unnormalized density weights
  vec gr_tmp;
  if (deriv)
  {
    gr = zeros(n_shortpar);  // gradient
    he = zeros(n_shortpar, n_shortpar); // hessian
  }
  for (uword q = 0; q < quad_y->n_quad; ++q)  // compute densities
  {
    udns(q) = trunc_exp( basis_exp(gr_tmp, quad_y->node(q, 0), x, deriv) ) * 
      quad_y->weight(q);
    if (deriv) 
    {
      vec gr1 = gr_tmp * udns(q);
      gr += gr1;
      he += gr1 * gr_tmp.t();
    }
  }
  double log_norm_const = trunc_log( accu(udns) );  // log normalizing constant
  if (deriv)  // update derivatives
  {
    double norm_const_1 = trunc_exp(-log_norm_const);
    gr *= norm_const_1;
    he *= norm_const_1; 
    he -= gr * gr.t();
  }
  return log_norm_const;
}

/* cond_log_dns: conditional log density
 *
 * return: log density values (double, dim = y.n_elem x x.n_elem) */

mat Item::cond_log_dns(
  vec y,  // y values (double, dim = y.n_elem)
  mat x  // x values (double, dim = x.n_elem x dim)
  )
{
  // initialize
  vec gr;
  mat f(y.n_elem, x.n_rows), he;
  // loop over x and y values
  #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
  #endif
  for (uword j = 0; j < x.n_rows; ++j)
  {
    double log_nc = log_normalize(gr, he, x(j, dim), false);
    for (uword i = 0; i < y.n_elem; ++i)
      f(i, j) = basis_exp(gr, y(i), x(j, dim), false) - log_nc;
  }
  return f;
}

/* penalize: evaluate penalty and its derivatives
 *
 * returns: basis expansion value (double)
 *
 * updates: 
 * gr: gradient (if deriv = true, cumulative; vec, dim = n_shortpar)
 * he: Hessian (if deriv = true, cumulative; mat, 
 *     dim = n_shortpar x n_shortpar) */

double Item::penalize(
  vec& gr,  // gradient vector (vec&, dim = n_shortpar)
  mat& he,  // Hessian matrix (mat&, dim = n_shortpar x n_shortpar)
  bool deriv  // whether derivatives are evaluated (bool)  
  )
{
  mat spmat(shortpar.memptr(), n_shortpar_y, n_shortpar_x + 1, false);
  // penalty term
  // univariate part
  double p = dot( spmat.col(0), pen_y * spmat.col(0) );
  // bivariate part
  for (uword j = 1; j <= n_shortpar_x; ++j) 
    p += dot( spmat.col(j), pen_y * spmat.col(j) );
  for (uword i = 0; i < n_shortpar_y; ++i) 
  {
    p += dot(spmat.submat(i, 1, i, n_shortpar_x), 
      spmat.submat(i, 1, i, n_shortpar_x) * pen_x);
  }
  // derivatives
  if (deriv)
  {
    // gradient
    mat gmat(gr.memptr(), n_shortpar_y, n_shortpar_x + 1, false);
    // univariate part
    gmat.col(0) += pen_y * spmat.col(0);
    // bivariate part
    gmat.tail_cols(n_shortpar_x) +=
      pen_y * spmat.tail_cols(n_shortpar_x) + 
      spmat.tail_cols(n_shortpar_x) * pen_x;

    // Hessian
    // univariate part
    he.submat(0, 0, n_shortpar_y - 1, n_shortpar_y - 1) += pen_y;
    // bivariate part
    // column-wise penalty matrix
    for (uword j = 1; j <= n_shortpar_x; ++j)
    {
      uword start = j * n_shortpar_y;
      uword end = start + n_shortpar_y - 1;
      he.submat(start, start, end, end) += pen_y;
    }
    // row-wise penalty matrix
    for (uword i = 0; i < n_shortpar_y; ++i)
    {
      uvec indx = regspace<uvec>(n_shortpar_y + i, 
        n_shortpar_y, n_shortpar - 1);
      he.submat(indx, indx) += pen_x;
    }
  }
  return 0.5 * p;
}

/* mloglik: M-step minus log-likelihood with penalty
 *
 * returns: log-likelihood value (double)
 *
 * updates: 
 * grad: gradient (if deriv = true, destroy first; double, dim = n_shortpar)
 * hess: Hessian (if deriv = true, destroy first; double, 
 *       dim = n_shortpar x n_shortpar) */

void Item::mloglik(
  bool deriv  // whether derivatives evaluated (bool)
  )
{
  // initialization
  f = 0.0;
  if (deriv)
  {
    grad.zeros();
    hess.zeros();
  }

  #ifdef _OPENMP
    #pragma omp parallel
  #endif
  {
    // local variables
    double f_local = 0.0;
    vec gr_local, grad_local; 
    mat he_local, hess_local; 
    if (deriv)
    {
      grad_local = zeros(n_shortpar);
      hess_local = zeros(n_shortpar, n_shortpar);
    }

    // main loop over quadrature for x
    #ifdef _OPENMP
      #pragma omp for schedule(static)
    #endif
    for (uword p = 0; p < quad_x.n_quad; ++p)
    {
      double estep_wt_p = 0.0;
      // person loop
      for (uword i = 0; i < dat.n_elem; ++i)
      {
        // skip if missing
        if ( is_equal(dat(i), na) ) continue;
        estep_wt_p += estep_wt(p, i);  // accumulate estep weights at quadrature p
        // log-likelihood
        f_local -= estep_wt(p, i) * 
          basis_exp(gr_local, dat(i), quad_x.node(p, dim), deriv);
        if (deriv) grad_local -= estep_wt(p, i) * gr_local;  // gradient
      }

      // normalizing constant (outside of person loop)
      log_norm_const(p) = 
        log_normalize(gr_local, he_local, quad_x.node(p, dim), deriv);
      f_local += estep_wt_p * log_norm_const(p);
      if (deriv)
      {
        grad_local += estep_wt_p * gr_local;
        hess_local += estep_wt_p * he_local;
      }
    }

    // accumulate
    #ifdef _OPENMP
      #pragma omp critical
    #endif
    {
      f += f_local;
      if (deriv)
      {
        grad += grad_local;
        hess += hess_local;
      }
    }
  }

  // penalty
  pen_val = penalize(grad, hess, deriv);
  f += pen_val;
}
