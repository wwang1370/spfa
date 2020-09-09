/* Constructors and methods: Class _Bspline_
 *
 * Author: Yang Liu
 *
 * Last modified: 09/01/2020 */

#include "bspline.h"

/* basic constructor */

Bspline::Bspline(
  uword n_basis_,  // number of basis functions (int)
  uword order_  // order of the spline (int) 
  ) : 
  n_basis(n_basis_), 
  order(order_)
{
  knots = eq_spc_knots();
}

/* constructor with side condition applied */

Bspline::Bspline(
  uword n_basis_,  // number of basis functions (int)
  uword order_,  // order of the spline (int) 
  double x0_,  // value at which side condition is applied (double) 
  double lmbd_ // penalty weights (double)
  ) : 
  Bspline(n_basis_, order_)
{
  // side condition (transformed, YL 09/01/2020)
  x0 = x0_;
  null0 = null_side().t();
  null0 = solve(null0 * diff_mat(1).t(), null0);
  inplace_trans(null0);
  // penalty
  lmbd = lmbd_;
  pen = pen_mat();
}

/* constructor including lower and upper envelopes */

Bspline::Bspline(
  uword n_basis_,  // number of basis functions (int)
  uword order_,  // order of the spline (int) 
  double x0_,  // value at which side condition is applied (double) 
  double lmbd_,  // penalty weights (double)
  uword n_sub_  // number of sub-intervals (double)
  ) : 
  Bspline(n_basis_, order_, x0_, lmbd_)
{
  // envelope
  n_sub = n_sub_;
  uword n_sub_total = n_sub * (n_basis - order + 1);
  lower.set_size(n_sub_total, n_basis);
  upper.set_size(n_sub_total, n_basis);
  env();
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
      // 2nd term
      denom = knots(which + order) - knots(which + 1);
      if (std::abs(denom) >= DBL_EPS)  // avoid dividing by 0
      {
        b2 = eval(x, which + 1, order - 1);
        t2 = ( knots(which + order) - x ) * b2 / denom;
      }
      f = t1 + t2;
    }
  }
  return f;
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

/* eval: evaluate B-spline basis at multiple data points
 *
 * returns: basis matrix (mat, dim = n x n_basis) */

mat Bspline::eval(
  const mat& x  // data points (const mat&)
  ) 
{
  mat basis(x.n_elem, n_basis);
  for (uword j = 0; j < n_basis; ++j) 
  {
    for (uword i = 0; i < x.n_elem; ++i)
      basis(i, j) = eval(x(i), j, order);
  }
  return basis;
}

/* eval0: evaluate B-spline basis at a single data point, with side condition
 *
 * returns: basis vector (rowvec, dim = n_basis) */

rowvec Bspline::eval0(
  double x  // data point (double)
  ) 
{
  rowvec basis = eval(x);
  return basis * null0;
}

/* eval0: evaluate B-spline basis at multiple data points, with side condition
 *
 * returns: basis matrix (mat, dim = n x n_basis) */

mat Bspline::eval0(
  const mat& x  // data points (const mat&)
  ) 
{
  mat basis = eval(x);
  return basis * null0;
}

/* eq_spc_knots: equally spaced knot sequence
 *
 * returns:
 * knots: knot sequence (vec, dim = n_basis + order) */

vec Bspline::eq_spc_knots() 
{
  double o_1 = (double) order - 1.0,            // order - 1
    d = 1.0 / ( (double) n_basis - o_1 );       // space between knots
  vec knots = linspace(- o_1 * d, 1.0 + o_1 * d, n_basis + order);
  return knots;
}


/* env: obtain upper and lower envelope for B-spline basis
 *
 * updates:
 * lower: lower envelope matrix (mat, dim = n_int x n_basis)
 * upper: upper envelope matrix (mat, dim = n_int x n_basis) */

void Bspline::env()
{
  vec grid = linspace(0.0, 1.0, lower.n_rows + 1);
  // loop over basis
  for (uword j = 0; j < lower.n_cols; ++j)
  {  
    uword prev;
    // first sub-interval for lower
    lower(0, j) = eval(grid(0), j, order);
    // middle sub-intervals
    for (uword i = 1; i < lower.n_rows; ++i)
    {
      lower(i, j) = eval(grid(i), j, order);
      prev = i - 1;
      upper(prev, j) = lower(i, j);
      if ( lower(prev, j) > upper(prev, j) )
        swap( lower(prev, j), upper(prev, j) );
    } 
    // last sub-interval for upper
    prev = lower.n_rows - 1;
    upper(prev, j) = eval(grid(lower.n_rows), j, order);
    if ( lower(prev, j) > upper(prev, j) )
      swap( lower(prev, j), upper(prev, j) );
  }
}

/* null_side: compute null space of side conditions
 *
 * returns: null space matrix (mat, dim = n_basis x (n_basis - 1)) */

mat Bspline::null_side()
{
  rowvec basis0 = eval(x0);   // evaluate basis function at x0
  mat null0 = null(basis0);   // null space of basis0
  return null0;
}

/* diff_mat: evaluate the difference matrix of a given order
 *
 * returns: differencing matrix (double, dim = (n_basis - order) x n_basis) */

mat Bspline::diff_mat(
  uword order  // order of difference (int)
  )
{
  mat diff(n_basis - order, n_basis); diff.zeros();
  if (order == 1)  // exit: order = 1
  {
    diff.diag(0).fill(1.0);
    diff.diag(1).fill(-1.0);
  }
  else  // otherwise, do recursion
  {
    uword o_1 = order - 1;
    // apply 1st-order diff matrix to columns of diff_mat(o_1)
    mat diff_ = diff_mat(o_1);
    for (uword j = 0; j < n_basis; ++j)
    {
      uword i = j < order ? 0 : j - order;
      while (i < j + 1 && i < diff.n_rows)
      {
        diff(i, j) = diff_(i, j) - diff_(i + 1, j);
        ++i;
      }
    }
  }
  return diff;
}

/* pen_mat: penalty matrix
 *
 * returns: penalty matrix ( mat, dim = (n_basis - 1) x (n_basis - 1) ) */

mat Bspline::pen_mat()
{
  mat d2null0 = diff_mat(2) * null0;   // transformed to penalize shortpar 
  mat pen = lmbd * d2null0.t() * d2null0;
  return pen;
}
