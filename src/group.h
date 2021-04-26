/* Headers: class Group
 *
 * Author: Yang Liu
 *
 * Note: Consider merging class Group with class Item later
 *
 * Last modified: 04/24/2021 */

#ifndef GROUP_H
#define GROUP_H

#include "basis.h"
#include "quad.h"

// class Group
class Group
{
  private:
    // variables to appear in the initialization list
    const uword n_dim = 2;  // # of dimensions (to be extended later)
    Bspline &basis_x;  // basis for x
    const mat &pen_x;  // penalty matrix for x
    const Quad &quad_x;  //  quadrature for x
    mat &estep_wt;  // E-step posterior weights for samples

    // intermediate variables and placeholders
    uword n_par;  // total number of parameters
    vec par;  // item parameters
    vec norm_const;  // normalizing constant for basis
    double f;  // M-step minus log-likelihood
    vec grad;  // gradient of M-step minus log-likelihood
    mat hess;  // Hessian of M-step minus log-likelihood
    vec dir;  // search direction
    double cond1;   // first-order condition 
    double pen_val;  // penalty value
    uvec activ;  // active constraints

    // methods
    double basis_exp( vec& gr, rowvec x, bool deriv);  // basis expansion
    double penalize(vec &gr, mat &he, bool deriv);   // evaluate penalty
    mat constr_mat();  // evaluate constraint matrix
    void mloglik(bool deriv);  // compute M-step minus log-likelihood
    void work_qp(vec &p, vec &mult); // solve the working QP
    void line_search();  // backtracking line search
    void search_dir();  // solution to the linearly constrained QP

  public:
    // methods
    void mstep(uword maxit, double tol);  // M-step optimization
    inline vec get_par() {return par;}  // retrieve par

    // constructor and destructor
    Group(Bspline &basis_x_, const mat &pen_x_, const Quad &quad_x_, 
      mat &estep_wt_);
    ~Group() {};

};

#endif
