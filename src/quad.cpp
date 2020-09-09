/* Constructors and methods: Class _Quad_
 *
 * Author: Yang Liu
 *
 * Last modified: 08/26/2020 */

#include "quad.h"

/* constructor */

Quad::Quad(
  uword n_quad_,  // number of quadrature nodes (int)
  double lower,  // lower limit (double)
  double upper  // upper limit (double)
  ) : 
  n_quad(n_quad_)
{
  gauss_legendre(lower, upper);
}

/* gauss_legendre: rescaled Gauss-Legendre quadrature
 *
 * returns: none
 *
 * updates:
 *   node: quadrature nodes (double, dim = 1)
 * weight: quadrature weights (double, dim = 1) */

void Quad::gauss_legendre(
  double lower,  // lower limit (double, dim = 1)
  double upper   // upper limit (double, dim = 1)
)
{
  vec k = linspace(1, n_quad - 1, n_quad - 1);
  mat jacobi = zeros(n_quad, n_quad);
  // beta diagonal from Golub & Welsch
  vec beta = k % arma::pow( (2 * k + 1) % (2 * k - 1) , -0.5 );
  // setting off-diagonal elements
  jacobi.diag(1) = beta;
  jacobi.diag(-1) = beta;
  // eigen decomposition
  vec eigen_val;
  mat eigen_vec;
  eig_sym(eigen_val, eigen_vec, jacobi);
  // determine nodeds and weights
  weight = 2 * vectorise( eigen_vec.row(0) % eigen_vec.row(0) );
  node = 0.5 * ( (upper - lower) * eigen_val + lower + upper );
  weight *= 0.5 * (upper - lower);
}
