/* Methods: EM algorithm
 *
 * Author: Yang Liu
 *
 * Last modified: 02/11/2021 */

#include "test.h"

/* search_dir0: Unconstrained Newton search direction
 *
 * update: search direction (vec, dim = n_shortpar) */ 

void Item::search_dir0()
{
  // solve for Newton direction if reciprocal condition number is not too small
  // (YL: 02/11/21)
  double kappa = rcond(hess);
  if (kappa > TOL_NEWT)
    dir = - solve(hess, grad, solve_opts::fast);
  // otherwise fall back to gradien ascent
  else
  {
    Rcout << "Warning: Hessian is ill-conditioned." << endl;
    dir = - grad;
  }
}

/* search_dir1: Solve the linearly constrained QP
 *
 * update: search direction (vec, dim = n_shortpar) */ 

void Item::search_dir1()
{
  vec p, mult;
  for (uword r = 0; r < MAX_QP; ++r)
  {
    // solve the working QP
    work_qp(p, mult);
    // if not moving
    if ( all(arma::abs(p) < DBL_EPS) )
    {
      if ( all(mult >= 0.0) ) break; // terminate the program
      uvec indx_activ = find(activ == 1);
      activ( indx_activ( mult.index_min() ) ) = 0;  // drop one constraint
    }
    else
    {
      uvec blk = find(activ == 0 && p < 0.0); // potential blocking constraints
      double alpha = 1.0;
      if (blk.n_elem > 0)
      {
        vec ratio = ( - shortpar.elem(blk) - dir.elem(blk) ) / p.elem(blk);
        uword indx_minr = ratio.index_min();
        if (ratio(indx_minr) < 1.0)
        {
          alpha = ratio(indx_minr);
          activ( blk(indx_minr) ) = 1;  // add one constraint
        }
      }
      dir += alpha * p;
    }
  }
}

/* work_qp: working equality constrained QP using the range-space approach
 *
 * update: 
 *    p: search sub-direction (vec, dim = n_shortpar)
 * mult: Lagrange multipler (vec, dim = n_activ) */ 

void Item::work_qp(
  vec& p,  // sub-direction (vec, dim = n_shortpar)
  vec& mult  // Lagrange multipler (vec, dim = n_activ)
  )
{
  p.set_size(n_shortpar); mult.set_size( size(activ) );
  uvec indx_activ = find(activ == 1), indx_inact = find(activ != 1);
  // sub-direction
  p.elem(indx_activ) = - dir(indx_activ) - shortpar.elem(indx_activ);
  p.elem(indx_inact) = - dir.elem(indx_inact) + 
    solve( hess.submat(indx_inact, indx_inact),
    hess.submat(indx_inact, indx_activ) * shortpar.elem(indx_activ) - 
    grad.elem(indx_inact) );
  // Lagrange multiplier
  mult = - hess.submat(indx_activ, indx_activ) * shortpar.elem(indx_activ) +
    hess.submat(indx_activ, indx_inact) * 
    ( dir.elem(indx_inact) + p.elem(indx_inact) ) + grad.elem(indx_activ);
}

/* line_search: backtracking line search
 *
 * update: 
 * shortpar: parameter vector (vec, dim = n_shortpar) */ 

void Item::line_search()
{
  double step_size = 1.0;
  double f0 = (this->f), inprod = - dot(dir, grad);
  vec shortpar0 = shortpar;
  for (uword b = 0; b < MAX_BKTRK; ++b)
  {
    shortpar = shortpar0 + step_size * dir;
    mloglik(false);
    if (f0 - (this->f) > ARMIJO2 * step_size * inprod) break;
    step_size *= ARMIJO1;
  }
}

/* mstep: item-level M-step 
 *
 * update: 
 *  shortpar: parameter vector (vec, dim = n_shortpar)
 * rough_pen: roughness penalty (double) */ 

void Item::mstep(
  uword maxit,  // maximum number of iterations (int)
  double tol  // convergence tolerance (double)
  )
{
  // initialize
  mloglik(true);
  // (proximal) Newton iterations
  for (uword r = 0 ; r < maxit; ++r)
  {
    (this->*search_dir_ptr)();  // search direction
    line_search();  // line search
    mloglik(true);  // update (this->f), grad, and hess
    if (arma::norm(grad) < tol) break;  // check convergence
  }
}

/* mstep: test-level M-step wrapper
 *
 * update: 
 * items[j].shortpar: parameter vector (vec, dim = n_shortpar) */ 

void Test::mstep()
{
  for (uword j = 0; j < n_item; ++j)
    items[j].mstep(maxit_mstep, tol_mstep);
}

/* estep: E-step
 *
 * update: 
 *  estep_wt: E-step weights (mat, dim = n_obsn x n_quad) 
 * (this->f): penalized marginal log-likelihood (double) */ 

void Test::estep()
{
  // initialize
  (this->f) = 0.0;
  estep_wt.zeros();
  vec gr;
  // person loop
  for (uword i = 0; i < n_obsn; ++i)
  {
    // item loop
    for (uword j = 0; j < n_item; ++j) 
    {
      // skip if missing
      if (items[j].dat(i) < 0 || items[j].dat(i) > 1) continue;
      // quadrature loop
      for (uword p = 0; p < quad.n_quad; ++p)
      {
        estep_wt(p, i) += 
          items[j].basis_exp(gr, items[j].dat(i), quad.node(p), false) - 
          items[j].log_norm_const(p);
      }
    }
    estep_wt.col(i) = arma::exp( estep_wt.col(i) ) % quad.weight;
    // marginal likelihood (add underflow protection, YL 02/08/21)
    double marg_lik = std::max(accu( estep_wt.col(i) ), MIN_ML);
    double marg_lik_1 = 1.0 / marg_lik;
    (this->f) += log(marg_lik);
    // normalize
    estep_wt.col(i) *= marg_lik_1;
  }
  // penalty
  for (uword j = 0; j < n_item; ++j) 
    (this->f) -= items[j].rough_pen;
}

/* start_val: starting values for item parameters
 *
 * update: 
 *  estep_wt: do not use in EM (vec, dim = n_shortpar)
 *  shortpar: parameter vector (vec, dim = n_shortpar)
 * rough_pen: roughness penalty (double) */ 

void Test::start_val()
{
  // re-configure estep_wt, maxit, and tol
  estep_wt.zeros();
  for (uword i = 0; i < n_obsn; ++i)
  {
    double x = 0.0;
    uword n = 0;
    for (uword j = 0; j < n_item; ++j)
    {
      if (items[j].dat(i) >= 0 && items[j].dat(i) <= 1)
      {
        x += items[j].dat(i);
        ++n;
      }
    }
    if (n > 0) x /= n;  // use observed mean as starting value for x
    estep_wt(arma::abs(quad.node - x).index_min(), i) = 1.0;
  }
  // run a single M-step
  for (uword j = 0; j < n_item; ++j)
  {
    Rcout << "Starting values: Item " << j << '\r';  // print info
    items[j].mstep(maxit_start, tol_start);
  }
  Rcout << endl;
}


/* em: Bock-Aitken EM algorithm
 *
 * update: 
 *  items[j].shortpar: item parameters (vec, dim = n_shortpar) 
 * items[j].rough_pen: roughness penalty (double)
 *          (this->f): penalized marginal log-likelihood (double) */ 

void Test::em()
{
  time_t begin, end;
  begin = clock();
  // main loop
  double f0 = DBL_MAX;
  for (iter = 0; iter < maxit_em; ++iter)
  {
    estep();  // estep
    Rcout << "EM iter " << iter << ": Penalized LL = " << 
      fixed << setprecision(4) << (this->f) << 
      '\r';  // print info
    if (std::abs( (this->f) - f0 ) < tol_em) break; // check convergence
    f0 = (this->f);
    mstep(); // M-step
  }
  Rcout << endl;   
  end = clock();
  time = double(end - begin) / double(CLOCKS_PER_SEC);  // record execution time
}
