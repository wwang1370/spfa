/* Constructor and Methods for Class Group
 *
 * Author: Yang Liu
 *
 * Last modified: 05/04/2021 */

#include "group.h"

/* constructor */
Group::Group(
  const vec &par_,  // parameters (vec), dim = n_par)
  Bspline &basis_x_,  // basis for x (Bspline)
  const mat &pen_x_,  // penalty matrix, x (mat)
  const Quad &quad_x_,  // quadrature for x (Quad)
  mat &estep_wt_  // E-step posterior weights (mat&, dim = n_obsn x n_quad)
  ) : 
  par(par_),
  basis_x(basis_x_), pen_x(pen_x_), 
  quad_x(quad_x_), estep_wt(estep_wt_)
{
  // initialize
  n_par = pow(basis_x.n_basis, n_dim);  // total number of parameters
  if (n_par != par.n_elem)
    throw runtime_error("n_par is not the same as the length of par.");

  // derivatives
  grad.set_size(n_par); 
  hess.set_size(n_par, n_par);
  dir.zeros(n_par);
  norm_const = basis_x.get_norm_const();
  cond1 = DBL_MAX;

  // active set
  activ.set_size(n_par);
  activ.fill(0);
}

/* basis_exp: evaluate basis expansion and its derivatives
 *
 * returns: basis expansion value (double)
 *
 * updates: 
 * gr: gradient (if deriv = true, destroy first; vec, dim = n_par) */

double Group::basis_exp(
  vec& gr,  // gradient vector (vec&, dim = n_shortpar)
  rowvec x,  // x values (rowvec)
  bool deriv  // whether derivatives evaluated (bool)
  )
{
  mat pmat(par.memptr(), basis_x.n_basis, basis_x.n_basis, false);
  // basis expansion
  rowvec b0 = basis_x.eval( x(0) ), b1 = basis_x.eval( x(1) );
  double f = dot(b0 * pmat, b1);
  // gradient
  if (deriv)
  {
    gr.set_size(n_par);
    gr = kron(b1, b0).t();
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

double Group::penalize(
  vec& gr,  // gradient vector (vec&, dim = n_shortpar)
  mat& he,  // Hessian matrix (mat&, dim = n_shortpar x n_shortpar)
  bool deriv  // whether derivatives are evaluated (bool)  
  )
{
  mat pmat(par.memptr(), basis_x.n_basis, basis_x.n_basis, false);
  double p = 0.0;
  // penalty term
  for (uword j = 0; j < basis_x.n_basis; ++j) 
    p += dot( pmat.col(j), pen_x * pmat.col(j) );  // columnwise
  for (uword i = 0; i < basis_x.n_basis; ++i) 
    p += dot(pmat.row(i), pmat.row(i) * pen_x);  // rowwise

  // derivatives
  if (deriv)
  {
    // gradient
    mat gmat(gr.memptr(), basis_x.n_basis, basis_x.n_basis, false);
    gmat += pen_x * pmat + pmat * pen_x;

    // Hessian
    // column-wise
    for (uword j = 0; j < basis_x.n_basis; ++j)
    {
      uword start = j * basis_x.n_basis;
      uword end = start + basis_x.n_basis - 1;
      he.submat(start, start, end, end) += pen_x;
    }
    // row-wise
    for (uword i = 0; i < basis_x.n_basis; ++i)
    {
      uvec indx = regspace<uvec>(i, basis_x.n_basis, n_par - 1);
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
 * grad: gradient (if deriv = true, destroy first; double, dim = npar)
 * hess: Hessian (if deriv = true, destroy first; double, 
 *       dim = npar x npar) */

void Group::mloglik(
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
    // define local variables
    double f_local = 0.0;
    vec gr_local, grad_local;
    mat hess_local;
    if (deriv)
    {
      grad_local = zeros(n_par);
      hess_local = zeros(n_par, n_par);
    }

    // main loop over quadrature for x
    #ifdef _OPENMP
      #pragma omp for schedule(static)
    #endif
    for (uword p = 0; p < quad_x.n_quad; ++p)
    {
      double loglik = trunc_log( 
        basis_exp(gr_local, quad_x.node.row(p), deriv) );
      double estep_wt_p = accu( estep_wt.row(p) );
      f_local -= estep_wt_p * loglik;  // update minus log-likelihood
      if (deriv)
      {
        double lik_1 = trunc_exp(-loglik);
        grad_local -= (estep_wt_p * lik_1) * gr_local;
        hess_local += (estep_wt_p * lik_1 * lik_1) * gr_local * gr_local.t();
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

/* constr_mat: constraint matrix
 * 
 * return: constraint matrix (mat, dim = n_constr x n_par) */

mat Group::constr_mat()
{
  // initialize
  uword n_activ = accu(activ),  // number of active inequality constraints
    n_eq = 2 * basis_x.n_basis - 1;  // number of equality constraints
  mat cmat = zeros(n_eq + n_activ, n_par);
  
  // equality constraints
  cmat.head_rows(basis_x.n_basis)  = 
    kron( eye(basis_x.n_basis, basis_x.n_basis), norm_const.t() );  // column-wise
  mat tmp = kron( norm_const.t(), eye(basis_x.n_basis, basis_x.n_basis) );
  cmat.rows(basis_x.n_basis, 2 * basis_x.n_basis - 2) = tmp.rows(0, 
    basis_x.n_basis - 2);  // row-wise

  // active inequality constraints
  uvec which_activ = find(activ == 1);
  for (uword k = 0; k < n_activ; ++k)
    cmat( n_eq + k, which_activ(k) ) = 1.0;  // active constraints
  return cmat;
}
