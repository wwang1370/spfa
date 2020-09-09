/* Constructors and basic methods: Class _test_
 *
 * Author: Yang Liu
 *
 * Last modified: 08/29/2020 */

#include "test.h"

/* constructor */

Test::Test(
  const mat& start,  // starting value matrix (mat, dim = n_shortpar x n_item)
  mat& dat_,  // data matrix (mat, dim = n_obsn x n_item)
  const uvec& type,   // item type (vec, dim = n_item)
  uword n_basis,  // number of basis functions (int)  
  double lmbd,  // penalty weight (double)
  uword n_quad,  // number of quadrature points (int) 
  uword maxit_em_,  // maximum number of EM iterations (int)
  uword maxit_mstep_,  // maximum number of M-step iterations (int)
  uword maxit_start_,  // maximum number of starting value iterations (int)
  double tol_em_,  // convergence tolerance for EM iterations (double)
  double tol_mstep_,  // convergence tolerance for M-step (double)
  double tol_start_  // convergence tolerance for starting values (double)
  ) : 
  n_obsn(dat_.n_rows), 
  n_item(dat_.n_cols),
  maxit_em(maxit_em_),
  maxit_mstep(maxit_mstep_),
  maxit_start(maxit_start_),
  tol_em(tol_em_),
  tol_mstep(tol_mstep_),
  tol_start(tol_start_)
{
  // initialize
  bspl = Bspline(n_basis, 4, 0.5, lmbd);
  quad = Quad(n_quad, 0.0, 1.0);
  estep_wt = mat(n_quad, n_obsn);
  // items
  for (uword j = 0; j < n_item; ++j)
  {
    Item item(start.col(j), dat_.colptr(j), n_obsn, type(j), bspl, quad, estep_wt,
      maxit_start);
    items.push_back(item);
  }
}

/* output: output routine
 *
 * return: list of objects to forward to R
 * 
 * args: none */

List Test::output()
{  
  // copy item parameters
  mat par(items[0].n_par, n_item);
  mat shortpar(items[0].n_shortpar, n_item);
  for (uword j = 0; j < n_item; ++j)
  {
    items[j].par.set_size(items[j].n_par);
    items[j].extend_par();  // shortpar -> par
    par.col(j) = items[j].par;
    shortpar.col(j) = items[j].shortpar;
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

