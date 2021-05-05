/* Headers: class Test
 *
 * Author: Yang Liu
 *
 * Last modified: 04/30/2021 */

#ifndef TEST_H
#define TEST_H

#include "group.h"
#include "item.h"

// constants
static const uword MAX_BKTRK = 10;   // max number of backtracking iterations
static const double ARMIJO1 = 0.5;   // first constant in line search
static const double ARMIJO2 = 0.0001;   // second constant in line search
static const uword MAX_QP = 20;   // max number of QP iterations
static const double TOL_SV = 1e-8;   // tolerance for rank determination
static const double TOL_SQP = 1e-8;   // tolerance for sequential QP

// class _Test_
class Test
{
  private:
    // variables
    const mat &dat;  // reference to data
    double na;  // missing code
    uword n_obsn, n_item;  // number of obsns and items
    uword n_dim;  // latent dimensionality
    uword maxit_em, maxit_mstep, 
      maxit_start;  // maximum number of EM, M-step, and starting value iterations
    double tol_em, tol_mstep,
      tol_start;  // tolerance for EM, M-step, and starting value iterations
    Bspline basis_x;  // common basis for both x and y
    GaussLegendre quad_x;  // quadrature for x
    //mat trans_x, pen_x_item, pen_x_group;  // common transformation and penalty matrices
    vector<Item> items;  // vector of items
    Group *group;  // group
    mat estep_wt;  // E-step weights
    uword iter;  // iteration counter
    double time;  // time consumed 
    int n_thrd;  // number of threads
    bool update_group;  // whether to update group parameters

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
    double risk(uvec it);  // L2 risk function

    // constructor and distructor
    Test(const mat &dat_, double na_, const uvec &item_type, 
      const Rcpp::List &start, const Rcpp::List &pos, 
      uword n_basis, vec lmbd, 
      uword n_quad, const uvec &dim, bool update_group_,
      uword maxit_em_, uword maxit_m_, uword maxit_start_, 
      double tol_em_, double tol_m_, double tol_start_,
      int n_thrd);
    ~Test();

};

#endif
