/* Headers: class _item_
 *
 * Author: Yang Liu
 *
 * Last modified: 08/27/2020 */

#ifndef ITEM_H
#define ITEM_H

#include "bspline.h"
#include "quad.h"

// class _Item_
class Item
{
  private:
    // variables
    uword n_par, n_shortpar;  // length of par and shortpars
    vec dat;  // data vector
    const uword type;  // item type (0 = no constraint, 1 = monotone)
    double f;  // M-step minus log-likelihood
    vec grad;  // gradient of M-step minus log-likelihood
    mat hess;  // Hessian of M-step minus log-likelihood
    vec dir;  // search direction
    double rough_pen;  // roughness penalty value
    vec log_norm_const;  // log normalizing constants
    uvec activ;  // active constraints

    // references
    Bspline& bspl;   // reference to B-spline bases
    Quad& quad;  // reference to quadrature
    mat& estep_wt;  // E-step posterior weights for samples

    // functions
    double basis_exp(vec& gr, double y, double x, bool deriv); // basis expansion
    double log_normalize(vec& gr, mat& he, double x, 
      bool deriv);   // evaluate log normalizing constant
    double pen(vec& gr, mat& he, bool deriv);   // evaluate penalty
    void mloglik(bool deriv);  // compute M-step minus log-likelihood
    void work_qp(vec& p, vec& mult); // solve the working QP
    void line_search();  // backtracking line search
    void mstep(uword maxit, double tol);  // M-step optimization
    
    // function pointers and target functions
    void (Item::*search_dir_ptr)();  // search direction 
    void search_dir0();  // solution to the Newton equation
    void search_dir1();  // solution to the linearly constrained QP

  public:
    // constructors
    // Note: Always get segfault when trying to pass a column of the data
    // matrix as a constant reference. As a workaround, pass dat_ptr as well as
    // n_obsn, and use the advanced constructor of vec in Item to avoid data
    // copy (YL 08/31/2020)
    Item(const vec& shortpar_, double* dat_ptr, uword n_obsn, uword type_, 
      Bspline& bspl_, Quad& quad_, mat& estep_wt_, uword maxit_start);

    // function
    void extend_par();  // shortpar -> par
    void reduce_par();  // par -> shortpar
    mat cond_log_dns(vec y, vec x);  // evaluate conditional log density

    // variables
    vec shortpar, par;  // item parameters

  friend class Test;
};

#endif
