/* Constructors and methods: Class Basis
 *
 * Author: Yang Liu
 *
 * Last modified: 04/13/2021 */

#include "basis.h"

/******************************************
 *          Derived Class Bspline
 *******************************************/

/* basic constructor */

Bspline::Bspline(
  uword n_basis_,  // number of basis functions (int)
  uword order_,  // order of the spline (int) 
  double lwr_,  // lower bound (double)
  double upr_  // upper bound (double)
  ) : Basis(n_basis_, lwr_, upr_)
{
  order = order_;  // order of B-spline
  knots = eq_spc_knots();  // equally-spaced knots
}

/* eval: evaluate the _which_th b-spline basis at a single datum
 *
 * the code is based on samiran sinha's implementation of de boor's recursive
 * algorithm, which can be found at:
 * 
 * http://www.stat.tamu.edu/~sinha/research/note1.pdf
 *
 * returns: basis function value (double) */

double Bspline::eval(
  double x,  // data value (double or vec)
  uword which,  // which basis (int) 
  uword order  // order of spline (int)
  )
{
  double denom, b1, b2, t1 = 0.0, t2 = 0.0, f = 0.0;
  // order = 1, piecewise linear function
  if (order == 1)
  {
    // localize to the support
    if (x >= knots(which) && x < knots(which + 1) )
      f = 1.0;
  }
  // otherwise, do recursion
  else
  {
    // localize to the support
    if ( x >= knots(which) && x <= knots(which + order) )
    {
      // 1st term
      denom = knots(which + order - 1) - knots(which);
      if (std::abs(denom) >= DBL_EPS)   // avoid dividing by 0
      {
        b1 = eval(x, which, order - 1);
        t1 = ( x - knots(which) ) * b1 / denom;
      }
      else t1 = 0.0;
      // 2nd term
      denom = knots(which + order) - knots(which + 1);
      if (std::abs(denom) >= DBL_EPS)  // avoid dividing by 0
      {
        b2 = eval(x, which + 1, order - 1);
        t2 = ( knots(which + order) - x ) * b2 / denom;
      }
      else t2 = 0.0;
      f = t1 + t2;
    }
  }
  return f;
}

/* eq_spc_knots: equally spaced knot sequence
 *
 * returns:
 * knots: knot sequence (vec, dim = n_basis + order) */

vec Bspline::eq_spc_knots() 
{
  double o_1 = (double) order - 1.0;  // order - 1
  double r = upr - lwr;  // range
  double d = r / ( (double) n_basis - o_1 );  // space between knots
  vec knots = linspace(lwr - o_1 * d, upr + o_1 * d, n_basis + order);
  return knots;
}

/* eval: evaluate B-spline basis at a single data point
 *
 * returns: basis vector (rowvec, dim = n_basis) */

rowvec Bspline::eval(
  double x  // data point (double)
  ) 
{
  rowvec basis(n_basis);
  for (uword j = 0; j < n_basis; ++j)
    basis(j) = eval(x, j, order);
  return basis;
}

/* get_norm_const: compute all normalizing constants
 *
 * returns: normalizing constants (vec, dim = n_basis) */

vec Bspline::get_norm_const()
{
  Bspline bs1(n_basis + 1, order + 1, lwr, upr);  // new basis
  vec delta = trans( bs1.eval(lwr) - bs1.eval(upr) );
  vec nc = cumsum(delta) * (upr - lwr) / ( (double) (n_basis - order + 1) );
  return nc.head(n_basis);
}

/******************************************
 *          Derived Class Iden
 *******************************************/

/* basic constructor */

Iden::Iden(
  uword n_basis_  // number of basis functions (int)
  ) : Basis(n_basis_) 
{
  lwr = 0.0;
  upr = (double) (n_basis - 1);
}

/* eval: evaluate identity basis at a single data point
 *
 * returns: basis vector (rowvec, dim = n_basis) */

rowvec Iden::eval(
  double x  // data point (double)
  ) 
{
  rowvec basis(n_basis); basis.zeros();
  uword indx( round(x) );
  if (indx >= lwr && indx <= upr)
    basis(indx) = 1.0;
  return basis;
}
