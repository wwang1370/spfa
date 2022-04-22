/* Constructors and methods: Class _Quad_
 *
 * Author: Yang Liu
 *
 * Last modified: 04/21/2022 */

#include "quad.h"

/******************************************
 *            Base Class Quad
 *******************************************/

/* basic constructor */

Quad::Quad(
  uword n_quad1_,  // number of quadrature points per dimension (uword)
  uword n_dim_  // dimension (uword);
  ) : n_quad1(n_quad1_), n_dim(n_dim_)
{
  n_quad = pow_uword(n_quad1, n_dim);
  node.set_size(n_quad, n_dim);
  node.set_size(n_quad, n_dim);
}

/* constructor with lower and upper bounds */

Quad::Quad(
  uword n_quad1_,  // number of quadrature points per dimension (uword)
  uword n_dim_,  // dimension (uword);
  double lwr_,  // lower limit (double)
  double upr_  // upper limit (double)
  ) : Quad(n_quad1_, n_dim_)
{
  lwr = lwr_;
  upr = upr_;
}

/******************************************
 *       Derived Class GaussLegendre
 *******************************************/

/* constructor */

GaussLegendre::GaussLegendre(
  uword n_quad1_,  // number of quadrature nodes per dimension (int)
  uword n_dim_,  // dimension (uword);
  double lwr_,  // lower limit (double)
  double upr_  // upper limit (double)
  ) : Quad(n_quad1_, n_dim_, lwr_, upr_)
{
  vec k = linspace(1, n_quad1 - 1, n_quad1 - 1);
  mat jacobi = zeros(n_quad1, n_quad1);
  // beta diagonal from Golub & Welsch
  vec beta = k % arma::pow( (2 * k + 1) % (2 * k - 1) , -0.5 );
  // setting off-diagonal elements
  jacobi.diag(1) = beta;
  jacobi.diag(-1) = beta;
  // eigen decomposition
  vec eigen_val;
  mat eigen_vec;
  eig_sym(eigen_val, eigen_vec, jacobi);
  // determine nodes and weights for one dimension
  vec weight1 = 2 * vectorise( eigen_vec.row(0) % eigen_vec.row(0) );
  vec node1 = 0.5 * ( (upr - lwr) * eigen_val + lwr + upr );
  weight1 *= 0.5 * (upr - lwr);
  // expand to grid
  node = expand_grid(node1, n_dim);
  weight = prod( expand_grid(weight1, n_dim), 1 );
}

/******************************************
 *       Derived Class Const
 *******************************************/

/* constructor */

Const::Const(
  uword n_quad1_  // number of quadrature nodes per dimension (int)
  ) : Quad(n_quad1_)
{
  lwr = 0.0;  // lower bound = 0
  upr = (double) (n_quad1 - 1);  // upper bound = n_quad1 - 1
  weight = ones(n_quad1);
  node = regspace(lwr, upr);
}

