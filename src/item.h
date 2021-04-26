/* Headers: class Item
 *
 * Author: Yang Liu
 *
 * Last modified: 04/14/2021 */

#ifndef ITEM_H
#define ITEM_H

#include "basis.h"
#include "quad.h"

// class Item
class Item
{
  private:
    // variables to appear in the initialization list
    const vec dat;  // data vector
    double na;  // missing data code
    vec shortpar;  // item parameters
    const uvec pos;  // positivity constraints
    const uword dim;  // dimension indicator
    Basis &basis_x;  // basis for x
    const mat &trans_x;  // transformation matrix for x
    const mat &pen_x;  // penalty matrix for x
    const Quad &quad_x;  //  quadrature for x
    mat &estep_wt;  // E-step posterior weights for samples

    // variables to be manually initialized
    Basis *basis_y;  // basis for y
    mat trans_y;  // transformation matrix for y
    mat pen_y;  // penalty matrix for y
    Quad *quad_y;  // quadrature for y

    // intermediate variables and placeholders
    uword n_par, n_par_y, n_par_x;  // lengths related to par
    uword n_shortpar, n_shortpar_y, n_shortpar_x;  // lengths related to shortpar
    vec par;  // long parameters
    double f;  // M-step minus log-likelihood
    vec grad;  // gradient of M-step minus log-likelihood
    mat hess;  // Hessian of M-step minus log-likelihood
    vec dir;  // search direction
    double cond1;   // first-order condition 
    double pen_val;  // penalty value
    vec log_norm_const;  // log normalizing constants
    uvec activ;  // active constraints

    // methods
    double log_normalize(vec &gr, mat &he, double x, 
      bool deriv);   // evaluate log normalizing constant
    double penalize(vec &gr, mat &he, bool deriv);   // evaluate penalty
    void mloglik(bool deriv);  // compute M-step minus log-likelihood
    void work_qp(vec &p, vec &mult); // solve the working QP
    void line_search();  // backtracking line search
    void (Item::*search_dir_ptr)();  // search direction 
    void search_dir0();  // solution to the Newton equation
    void search_dir1();  // solution to the linearly constrained QP

  public:
    // methods
    inline double get_na() {return na;}  // retrieve na
    inline vec get_shortpar() {return shortpar;}  // retrieve shortpar
    inline vec get_par() {return par;}  // retrieve par
    inline uword get_dim() {return dim;}  // retrieve dim
    inline double get_log_norm_const(uword which) {
      return log_norm_const(which);}  // retrieve log-normalizing constant
    inline double get_pen_val() {return pen_val;}  // retrieve par
    inline void set_shortpar(vec shortpar_) {shortpar = shortpar_;}  // set shortpar
    inline void set_par(vec par_) {par = par_;}  // set par
    void extend_par();  // shortpar -> par
    void reduce_par();  // par -> shortpar
    double basis_exp(vec &gr, double y, double x, bool deriv); // basis expansion
    void mstep(uword maxit, double tol);  // M-step optimization
    mat cond_log_dns(vec y, mat x);  // evaluate conditional log density

    // constructor and destructor
    Item(const vec &dat_, double na_, 
      const uword type, const vec &shortpar_, const uvec &pos_, 
      const uword dim_,
      Basis &basis_x_, const mat &trans_x_, const mat &pen_x_,
      const Quad &quad_x, 
      mat &estep_wt_);
    Item(const Item &item);  // copy
    ~Item();

};

#endif
