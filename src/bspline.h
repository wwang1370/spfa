/* Headers: class _Bspline_
 *
 * Author: Yang Liu
 *
 * Last modified: 05/14/2020 */

#ifndef BSPLINE_H
#define BSPLINE_H

#include "RcppArmadillo.h"
using namespace std;
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// constants
static const double DBL_EPS = 1e-64;   // tolerance for comparison between floats

// class _Bspline_
class Bspline
{
  private:
    // variables
    uword n_basis;  // number of basis functions
    uword order;  // order of spline
    vec knots;  // knot sequence
    uword n_sub;  // number of sub-intervals between knots
    mat lower, upper;  // lower and upper envelope
    double x0;  // value at which side condition applies
    mat null0;  // null space of side condition at x0
    double lmbd;  // penalty weight
    mat pen;  // penalty matrix
    mat dnull0;  // diff_mat(1) * null0
    mat c_mat;  // left-hand side matrix of the linear constraint

    // functions
    vec eq_spc_knots();  // equally-spaced knots
    mat null_side();  // null space of side condition at x = x0
    mat diff_mat(uword order);  // difference matrix
    mat pen_mat();  // penalty matrix
    double eval(double x, uword which, 
      uword order);  // single datum, single basis, recursive
    rowvec eval(double x);  // single datum
    rowvec eval0(double x);  // single datum, with side condition
    void env();  // obtain lower and upper envelopes

  public:
    // constructors
    Bspline(){}  // default
    Bspline(uword n_basis_, uword order_);  // minimum
    Bspline(uword n_basis_, uword order_, 
      double x0, double lmbd);   // with side condition and penalty
    Bspline(uword n_basis_, uword order_, 
      double x0, double lmbd, uword n_sub);  // with upper and lower envelope

    // functions
    mat eval(const mat& x);  // multiple data
    mat eval0(const mat& x);  // multiple data, with side condition

  friend class Item;
  friend class Test;
};

#endif
