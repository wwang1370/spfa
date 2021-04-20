/* Headers: class Test
 *
 * Author: Yang Liu
 *
 * Last modified: 04/17/2021 */

#ifndef TEST_H
#define TEST_H

#include "item.h"
#include "util.h"

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
    const mat &dat;  // reference to data
    double na;  // missing code
    uword n_obsn, n_item;  // number of obsns and items
    uword maxit_em, maxit_mstep, 
      maxit_start;  // maximum number of EM, M-step, and starting value iterations
    double tol_em, tol_mstep,
      tol_start;  // tolerance for EM, M-step, and starting value iterations
    Bspline basis_x;  // common basis for both x and y
    GaussLegendre quad_x;  // quadrature for x
    mat trans_x, pen_x;  // common transformation and penalty matrices
    vector<Item> items;  // vector of items
    mat estep_wt;  // E-step weights
    uword iter;  // iteration counter
    double time;  // time consumed 
    int n_thrd;  // number of threads

    // methods
    void init_estep_wt(const uvec& dim);  // initialize E-step weights

  public:
    // variables
    double f;  // penalized marginal log-likelihood

    // methods
    void estep();  // E-step
    void mstep();  // M-step wrapper
    void em();   // wrapper for EM algorithm
    mat score(uword mode);  // compute scores
    List output();  // generate output
    vec marg_lik(mat y, uvec it);  // marginal likelihood
    double risk(uvec it);  // L2 risk function

    // constructor and distructor
    Test(const mat &dat_, double na_, const uvec &item_type, 
      const Rcpp::List &start, const Rcpp::List &pos, 
      uword n_basis, double lmbd, 
      uword n_quad, const uvec &dim,
      uword maxit_em_, uword maxit_m_, uword maxit_start_, 
      double tol_em_, double tol_m_, double tol_start_,
      int n_thrd);
    ~Test() {};

};

#endif
