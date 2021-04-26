/* Constructor and Public Methods for Class Test
 *
 * Author: Yang Liu
 *
 * Last modified: 04/19/2021 */

#include "test.h"

/* constructor */

Test::Test(
  const mat& dat_,  // data matrix (mat, dim = n_obsn x n_item)
  double na_,  // missing code (double)
  const uvec& item_type,  // item types (uvec, dim = n_item)
  const Rcpp::List &start,  // starting values (List)
  const Rcpp::List &pos,   // positive constraints (List)
  uword n_basis,  // number of basis functions (int)
  double lmbd,  // penalty weight (double)
  uword n_quad,  // number of quadrature points (int)
  const uvec &dim,  // dimension indicator (uvec, dim = )
  uword maxit_em_,  // maximum number of EM iterations (int)
  uword maxit_mstep_,  // maximum number of M-step iterations (int)
  uword maxit_start_,  // maximum number of starting value iterations (int)
  double tol_em_,  // convergence tolerance for EM iterations (double)
  double tol_mstep_,  // convergence tolerance for M-step (double)
  double tol_start_,  // convergence tolerance for starting values (double)
  int n_thrd_  // number of threads
  ) : 
  dat(dat_), na(na_), n_obsn(dat_.n_rows), n_item(dat_.n_cols),
  maxit_em(maxit_em_), maxit_mstep(maxit_mstep_), maxit_start(maxit_start_),
  tol_em(tol_em_), tol_mstep(tol_mstep_), tol_start(tol_start_),
  basis_x(n_basis, 4, 0.0, 1.0), 
  quad_x(n_quad, dim.max() + 1, 0.0, 1.0)
{
  // parallel computing
  n_thrd = min( n_thrd_, omp_get_max_threads() );
  #ifdef _OPENMP
    omp_set_num_threads(n_thrd);
  #endif

  // basis and transformation matrices
  rowvec b0 = basis_x.eval(0.5);  // side value
  trans_x = null(b0).t();
  trans_x = - solve(trans_x * diff_mat(basis_x.n_basis, 1).t(), trans_x);
  inplace_trans(trans_x);  // transformation matrix
  mat d2tr = diff_mat(basis_x.n_basis, 2) * trans_x;
  pen_x = lmbd * d2tr.t() * d2tr;  // penalty matrix

  // initialize items
  init_estep_wt(dim);  // E-step weights
  for (uword j = 0; j < n_item; ++j)
  {
    Rcout << "Starting values: Item " << j << '\r';  // print info
    vec start_j = start[j];
    uvec pos_j = pos[j];
    items.emplace_back(dat.col(j), na,
      item_type(j), start_j, pos_j, dim(j),
      basis_x, trans_x, pen_x, quad_x, 
      estep_wt);
    items[j].mstep(maxit_start, tol_mstep);  // run M-step once to get starting values
  }
  Rcout << endl;
}

/* init_estep_wt: initiate E-step weights
 *
 * update: 
 * estep_wt: do not use in EM (vec, dim = n_shortpar)
 * shortpar: parameter vector (vec, dim = n_shortpar)
 *  pen_val: roughness penalty (double) */ 

void Test::init_estep_wt(
    const uvec &dim  // dimension indicator (uvec)
  )
{
  estep_wt = zeros(quad_x.n_quad, n_obsn);
  #ifdef _OPENMP
    #pragma omp parallel
  #endif
  {
    #ifdef _OPENMP
      #pragma omp for schedule(static)
    #endif
    for (uword i = 0; i < n_obsn; ++i)
    {
      rowvec x = zeros(1, quad_x.n_dim);
      uword n = 0;
      for (uword j = 0; j < n_item; ++j)
      {
        if ( is_equal(dat(i, j), na) ) continue;
        else
        {
          x( dim(j) ) += dat(i, j);
          ++n;
        }
      }
      if (n > 0) x /= n;  // use empirical mean as starting value for x
      mat delta = arma::abs(quad_x.node.each_row() - x);
      estep_wt(arma::sum(delta, 1).index_min(), i) = 1.0;
    }
  }
}

/* marg_lik: compute marginal likelihood
 *
 * return: marginal likelihood (double, dim = y.n_elem) */

vec Test::marg_lik(
  mat y,  // y values (double, n_cols = it.n_elem)
  uvec it  // item combination (int, dim = it.n_elem)
  )
{
  mat cdns = zeros(y.n_rows, quad_x.n_quad);
  for (uword k = 0; k < it.n_elem; ++k)  // accumulate log conditional density
    cdns += items[it(k)].cond_log_dns(y.col(k), quad_x.node);
  vec f = trunc_exp(cdns) * quad_x.weight;
  //cout << this->f << endl;
  return f;
}

///* score: compute EAP score after transformation to normal scale
// *
// * return: EAP score and posterior variance (double, dim = n_obsn x 2) */
//
//mat Test::score(
//  uword mode  // mode of scores (uint): 0 = uniform, 1 = normal
//  )
//{  
//  vec x(quad_x.n_quad);
//  // transform x by qnorm if mode = 1
//  if (mode == 1)
//  {
//    for (uword q = 0; q < quad_x.n_quad; ++q)
//      x(q) = R::qnorm(quad_x.node(q, 0), 0.0, 1.0, 1, 0);
//  }
//  else if (mode == 0) x = quad_x.node.col(0);  // no need to transform for mode = 0
//  mat ret(n_obsn, 2);
//  ret.col(0) = estep_wt.t() * x;   // mean
//  ret.col(1) = estep_wt.t() * (x % x) - ret.col(0) % ret.col(0); // variance
//  return ret;
//}
//
/* output: output routine
 *
 * return: list of objects to forward to R */

List Test::output()
{  
  List par, shortpar;
  // compute long parameters
  for (uword j = 0; j < n_item; ++j)
  {
    shortpar.push_back( arma2r( items[j].get_shortpar() ) );
    items[j].extend_par();  // shortpar -> par
    par.push_back( arma2r( items[j].get_par() ) );
  }
  // return list
  List ret = List::create(
    Named("par") = par,
    Named("shortpar") = shortpar,
    Named("f") = (this->f),
    Named("niter") = iter,
    Named("time") = time
  );
  return ret;
}
