/* Headers: class _Test_
 *
 * Author: Yang Liu
 *
 * Last modified: 09/29/2020 */

#ifndef TEST_H
#define TEST_H

#include "item.h"

// constants
static const uword MAX_BKTRK = 10;   // max number of backtracking iterations
static const double ARMIJO1 = 0.5;   // first constant in line search
static const double ARMIJO2 = 0.0001;   // second constant in line search
static const uword MAX_QP = 20;   // max number of QP iterations

// class _Test_
class Test
{
  private:
    // variables
    uword n_obsn, n_item;  // number of obsns and items
    uword maxit_em, maxit_mstep, 
      maxit_start;  // maximum number of EM, M-step, and starting value iterations
    double tol_em, tol_mstep,
      tol_start;  // tolerance for EM, M-step, and starting value iterations
    Bspline bspl;  // B-spline bases
    Quad quad;  // quadrature
    vector<Item> items;  // vector of items
    mat estep_wt;  // E-step weights
    uword iter;  // iteration counter
    double time;  // time consumed 

    // functions
    void mstep();  // M-step wrapper

  public:
    // constructors
    Test(const mat& start, mat& dat_, const uvec& type, 
      uword n_basis, double lmbd, uword n_quad, 
      uword maxit_em_, uword maxit_m_, uword maxit_start_, 
      double tol_em_, double tol_m_, double tol_start_);

    // variables
    double f;  // penalized marginal log-likelihood

    // function
    void start_val();  // obtain starting values
    void estep();  // E-step
    void em();   // wrapper for EM algorithm
    vec score();  // compute scores
    List output();  // generate output
};

#endif
