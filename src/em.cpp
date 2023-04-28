/* Methods: EM algorithm
 *
 * Author: Yang Liu
 *
 * Last modified: 08/23/2021 */

#include "test.h"

/* search_dir0: Unconstrained Newton search direction
 *
 * update: search direction (vec, dim = n_shortpar) */ 

void Item::search_dir0()
{
  // check condition number of hess (YL 08/23/21)
  double kappa = cond(hess);
  if (kappa <= 1e+8)
    dir = - solve(hess, grad, solve_opts::fast + solve_opts::likely_sympd);
  else
  {
    Rcout << "WARNING: Ill-conditioned Hessian; condition # = " << kappa << endl;
    dir = - grad;
  }
  cond1 = arma::norm(grad);
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
    if ( p.is_zero(TOL_SQP) )
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
  p.set_size(n_shortpar);
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
  // grad of Lagrangian
  double norm1 = arma::norm(grad.elem(indx_activ) - mult),
    norm2 = arma::norm( grad.elem(indx_inact) );
  cond1 = std::sqrt(norm1 * norm1 + norm2 * norm2);
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
 * shortpar: parameter vector (vec, dim = n_shortpar)
 *  pen_val: roughness penalty (double) */ 

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
    if (cond1 < tol) break;  // check convergence
  }
}

/* search_dir: Solve the linearly constrained QP (Method for Group)
 *
 * update: search direction (vec, dim = n_par) */ 

void Group::search_dir()
{
  dir.zeros();
  vec p, mult;
  for (uword r = 0; r < MAX_QP; ++r)
  {
    // solve the working QP
    work_qp(p, mult);
    // if not moving
    if ( p.is_zero(TOL_SQP) )
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
        vec ratio = ( - par.elem(blk) - dir.elem(blk) ) / p.elem(blk);
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

/* work_qp: working equality constrained QP (method for Group)
 *
 * update: 
 *    p: search sub-direction (vec, dim = n_par)
 * mult: Lagrange multipler (vec, dim = n_activ) */ 

void Group::work_qp(
  vec& p,  // sub-direction (vec, dim = n_par)
  vec& mult  // Lagrange multipler (vec, dim = n_activ)
  )
{
  uword n_activ = accu(activ);
  mat A = constr_mat();
  mat U, V; vec s;
  svd( U, s, V, A.t() );
  uword rk = accu(s >= TOL_SV); // rank
  mat Y = U.head_cols(rk), Z = U.tail_cols(n_par - rk);  // range and null spaces
  U.set_size(0, 0); V.set_size(0, 0); s.set_size(0);
  vec h = A * dir;
  if (n_activ > 0)
    h.tail(n_activ) += A.tail_rows(n_activ) * par;
  vec g = hess * dir + grad;
  mat AY = A * Y;
  vec p_range = Y * solve(AY, -h);
  mat ZtG = Z.t() * hess;
  vec p_null = Z * solve(ZtG * Z, - ZtG * p_range - Z.t() * g);
  p = p_range + p_null;
  vec mult_full = solve( AY.t(), Y.t() * (g + hess * p) );
  mult = mult_full.tail(n_activ);  // Lagrange multiplier
  cond1 = arma::norm(grad - A.t() * mult_full);
}

/* line_search: backtracking line search
 *
 * update: 
 * par: parameter vector (vec, dim = n_par) */ 

void Group::line_search()
{
  double step_size = 1.0;
  double f0 = (this->f), inprod = - dot(dir, grad);
  vec par0 = par;
  for (uword b = 0; b < MAX_BKTRK; ++b)
  {
    par = par0 + step_size * dir;
    mloglik(false);
    if (f0 - (this->f) > ARMIJO2 * step_size * inprod) break;
    step_size *= ARMIJO1;
  }
}

/* mstep: group-level M-step 
 *
 * update: 
 *      par: parameter vector (vec, dim = n_par)
 *  pen_val: roughness penalty (double) */ 

void Group::mstep(
  uword maxit,  // maximum number of iterations (int)
  double tol  // convergence tolerance (double)
  )
{
  // initialize
  mloglik(true);
  // (proximal) Newton iterations
  for (uword r = 0 ; r < maxit; ++r)
  {
    search_dir();  // search direction
    line_search();  // line search
    mloglik(true);  // update (this->f), grad, and hess
    if (cond1 < tol) break;  // check convergence
  }
}

/* mstep: test-level M-step wrapper
 *
 * update: 
 * items[j].shortpar: parameter vector (vec, dim = n_shortpar)
 *         group.par: parameter vector (vec, dim = n_par) */ 

void Test::mstep()
{
  for (uword j = 0; j < n_item; ++j)
    items[j].mstep(maxit_mstep, tol_mstep);
  if (update_group)
    group->mstep(maxit_mstep, tol_mstep);
}

/* estep: E-step
 *
 * update: 
 *  estep_wt: E-step weights (mat, dim = n_obsn x n_quad) 
 * (this->f): penalized marginal log-likelihood (double) */ 

void Test::estep()
{
  // initialize
  f = 0.0;
  estep_wt.zeros();
  vec gr;

  // append density
  vec wt = quad_x.weight;
  if (update_group)
  {
    #ifdef _OPENMP
      #pragma omp parallel for schedule(static)
    #endif
    for (uword p = 0; p < quad_x.n_quad; ++p)
      wt(p) *= group->basis_exp(gr, quad_x.node.row(p), false);
  }

  // person loop
  #ifdef _OPENMP
    #pragma omp parallel
  #endif
  { 
    // local variables
    double f_local = 0.0;

    // parallel loop
    #ifdef _OPENMP
      #pragma omp for schedule(static) 
    #endif
    for (uword i = 0; i < n_obsn; ++i)
    {
      // item loop
      for (uword j = 0; j < n_item; ++j) 
      {
        // skip if missing
        if ( is_equal( dat(i, j), items[j].get_na() ) ) continue;
        // quadrature loop
        for (uword p = 0; p < quad_x.n_quad; ++p)
        {
          estep_wt(p, i) += 
            items[j].basis_exp(gr, dat(i, j), 
            quad_x.node( p, items[j].get_dim() ), false) - 
            items[j].get_log_norm_const(p);
        }
      }

      //// append density
      //if (update_group)
      //{
      //  for (uword p = 0; p < quad_x.n_quad; ++p)
      //  {
      //    estep_wt(p, i) += 
      //      trunc_log( group->basis_exp(gr, quad_x.node.row(p), false) );
      //  }
      //}

      // append quadrature weights
      estep_wt.col(i) = trunc_exp( estep_wt.col(i) ) % wt;

      // marginal likelihood
      double marg_lik = accu( estep_wt.col(i) );
      double marg_lik_1 = 1.0 / marg_lik;
      f_local += trunc_log(marg_lik);
      // normalize
      estep_wt.col(i) *= marg_lik_1;
    }

    // accumulate local variables
    #ifdef _OPENMP
      #pragma omp critical
    #endif
    {
      f += f_local;
    }
  }

  // penalty
  for (uword j = 0; j < n_item; ++j) 
    f -= items[j].get_pen_val();
}


/* em: Bock-Aitken EM algorithm
 *
 * update: 
 * items[j].shortpar: item parameters (vec, dim = n_shortpar) 
 *  items[j].pen_val: roughness penalty (double)
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
