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
  mat spmat(shortpar.memptr(), bspl.n_basis - 1, bspl.n_basis, false);
  rowvec basis_y = bspl.eval0(y), basis_x = bspl.eval0(x);
  // basis expansion
  double f = dot(spmat.col(0), basis_y) + 
    dot(basis_y * spmat.tail_cols(bspl.n_basis - 1), basis_x);
  // gradient
  if (deriv)
  {
    gr.set_size(n_shortpar);
    mat gmat(gr.memptr(), bspl.n_basis - 1, bspl.n_basis, false);
    gmat.col(0) = trans(basis_y);
    gmat.tail_cols(bspl.n_basis - 1) = gmat.col(0) * basis_x;
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
  vec udns(quad.n_quad);  // unnormalized density weights
  vec gr_tmp;
  if (deriv)
  {
    gr = zeros(n_shortpar);  // gradient
    he = zeros(n_shortpar, n_shortpar); // hessian
  }
  for (uword q = 0; q < quad.n_quad; ++q)  // compute densities
  {
    udns(q) = exp( basis_exp(gr_tmp, quad.node(q), x, deriv) ) * quad.weight(q);
    if (deriv) 
    {
      vec gr1 = gr_tmp * udns(q);
      gr += gr1;
      he += gr1 * gr_tmp.t();
    }
  }
  double log_norm_const = log( accu(udns) );  // log normalizing constant
  if (deriv)  // update derivatives
  {
    double norm_const_1 = exp(-log_norm_const);
    gr *= norm_const_1;
    he *= norm_const_1; 
    he -= gr * gr.t();
  }
  return log_norm_const;
}

/* pen: evaluate penalty and its derivatives
 *
 * returns: basis expansion value (double)
 *
 * updates: 
 * gr: gradient (if deriv = true, cumulative; vec, dim = n_shortpar)
 * he: Hessian (if deriv = true, cumulative; mat, 
 *     dim = n_shortpar x n_shortpar) */

double Item::pen(
  vec& gr,  // gradient vector (vec&, dim = n_shortpar)
  mat& he,  // Hessian matrix (mat&, dim = n_shortpar x n_shortpar)
  bool deriv  // whether derivatives are evaluated (bool)  
  )
{
  mat spmat(shortpar.memptr(), bspl.n_basis - 1, bspl.n_basis, false);
  // penalty term
  // univariate part
  double p = dot( spmat.col(0), bspl.pen * spmat.col(0) );
  // bivariate part
  for (uword j = 1; j < bspl.n_basis; ++j) 
    p += dot( spmat.col(j), bspl.pen * spmat.col(j) );
  for (uword i = 0; i < bspl.n_basis - 1; ++i) 
  {
    p += dot(spmat.submat(i, 1, i, bspl.n_basis - 1), 
      spmat.submat(i, 1, i, bspl.n_basis - 1) * bspl.pen);
  }
  // derivatives
  if (deriv)
  {
    // gradient
    mat gmat(gr.memptr(), bspl.n_basis - 1, bspl.n_basis, false);
    // univariate part
    gmat.col(0) += bspl.pen * spmat.col(0);
    // bivariate part
    gmat.tail_cols(bspl.n_basis - 1) +=
      bspl.pen * spmat.tail_cols(bspl.n_basis - 1) + 
      spmat.tail_cols(bspl.n_basis - 1) * bspl.pen;

    // Hessian
    // univariate part
    he.submat(0, 0, bspl.n_basis - 2, bspl.n_basis - 2) += bspl.pen;
    // bivariate part
    // column-wise penalty matrix
    for (uword j = 1; j < bspl.n_basis; ++j)
    {
      uword start = j * (bspl.n_basis - 1);
      uword end = start + bspl.n_basis - 2;
      he.submat(start, start, end, end) += bspl.pen;
    }
    // row-wise penalty matrix
    for (uword i = 0; i < bspl.n_basis - 1; ++i)
    {
      uvec indx = regspace<uvec>(bspl.n_basis - 1 + i, 
        bspl.n_basis - 1, n_shortpar - 1);
      he.submat(indx, indx) += bspl.pen;
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
  (this->f) = 0.0;
  vec gr(n_shortpar);
  mat he(n_shortpar, n_shortpar);
  if (deriv)
  {
    grad.zeros();
    hess.zeros();
  }
  // main loop over quadrature for x
  for (uword p = 0; p < quad.n_quad; ++p)
  {
    // person loop
    double estep_wt_p = 0.0;
    for (uword i = 0; i < dat.n_elem; ++i)
    {
      // skip if missing
      if (dat(i) < 0.0 || dat(i) > 1.0) continue;
      estep_wt_p += estep_wt(p, i);
      // log-likelihood
      (this->f) -= estep_wt(p, i) * basis_exp(gr, dat(i), quad.node(p), deriv);
      if ( gr.has_nan() )
      {
        Rcout << "gr has nan" << endl;
        Rcout << shortpar.t() << endl;
      }
      if (deriv) grad -= estep_wt(p, i) * gr;  // gradient
    }
    // normalizing constant (outside of person loop)
    log_norm_const(p) = log_normalize(gr, he, quad.node(p), deriv);
    (this->f) += estep_wt_p * log_norm_const(p);
    if (deriv)
    {
      grad += estep_wt_p * gr;
      hess += estep_wt_p * he;
    }
  }
  // penalty
  rough_pen = pen(grad, hess, deriv);
  (this->f) += rough_pen;
}
