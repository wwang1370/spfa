/* Methods: Item models
 *
 * Author: Yang Liu
 *
 * Last modified: 08/28/2020 */

#include "item.h"

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
  vec x  // x values (double, dim = x.n_elem)
  )
{
  // initialize
  vec gr;
  mat f(y.n_elem, x.n_elem), he;
  // loop over x and y values
  for (uword j = 0; j < x.n_elem; ++j)
  {
    double log_nc = log_normalize(gr, he, x(j), false);
    for (uword i = 0; i < y.n_elem; ++i)
      f(i, j) = basis_exp(gr, y(i), x(j), false) - log_nc;
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
