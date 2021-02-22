/* Interfaces for internal classes
 *
 * Author: Yang Liu
 *
 * Last modified: 02/11/2021 */

#include "test.h"

/* risk: compute L2 risk function
 *
 * return: risk function (double) */

double Test::risk(
  uvec it  // item combination (int, dim = it.n_elem)
  )
{
  // 1st term
  double f = 0.0;
  // loop over the grid of quadrature points
  for (uword i = 0; i < pow(quad.n_quad, it.n_elem); ++i)
  {
    uvec qu = grid_loc(i, it.n_elem, quad.n_quad);  // grid location
    double mlik = as_scalar( marg_lik(quad.node(qu).t(), it) );
    f += mlik * mlik * arma::prod( quad.weight(qu) );
  }
  // 2nd term
  vec mlik = marg_lik(dat.cols(it), it);
  f -= 2.0 * arma::mean(mlik);
  return f;
}

/* marg_lik: compute marginal likelihood
 *
 * return: marginal likelihood (double, dim = y.n_elem) */

vec Test::marg_lik(
  mat y,  // y values (double, n_cols = it.n_elem)
  uvec it  // item combination (int, dim = it.n_elem)
  )
{
  mat cdns = zeros(y.n_rows, quad.n_quad);
  for (uword k = 0; k < it.n_elem; ++k)  // accumulate log conditional density
    cdns += items[it(k)].cond_log_dns(y.col(k), quad.node);
  vec f = trunc_exp(cdns) * quad.weight;
  return f;
}

/* score: compute EAP score after transformation to normal scale
 *
 * return: EAP score and posterior variance (double, dim = n_obsn x 2) */

mat Test::score(
  uword mode  // mode of scores (uint): 0 = uniform, 1 = normal
  )
{  
  vec x(quad.n_quad);
  // transform x by qnorm if mode = 1
  if (mode == 1)
  {
    for (uword q = 0; q < quad.n_quad; ++q)
      x(q) = R::qnorm(quad.node(q), 0.0, 1.0, 1, 0);
  }
  else if (mode == 0) x = quad.node;  // no need to transform for mode = 0
  mat ret(n_obsn, 2);
  ret.col(0) = estep_wt.t() * x;   // mean
  ret.col(1) = estep_wt.t() * (x % x) - ret.col(0) % ret.col(0); // variance
  return ret;
}

/* output: output routine
 *
 * return: list of objects to forward to R */

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

/* extend_par: shortpar -> par
 *
 * updates: par (double, dim = n_par) */

void Item::extend_par()
{
  mat spmat(shortpar.memptr(), bspl.n_basis - 1, bspl.n_basis, false),
    pmat(par.memptr(), bspl.n_basis, bspl.n_basis + 1, false);
  pmat.col(0) = bspl.null0 * spmat.col(0);
  pmat.tail_cols(bspl.n_basis) = bspl.null0 * 
    spmat.tail_cols(bspl.n_basis - 1) * bspl.null0.t();
}

/* reduce_par: par -> shortpar
 *
 * updates: shortpar (double, dim = n_shortpar) */

void Item::reduce_par()
{
  mat spmat(shortpar.memptr(), bspl.n_basis - 1, bspl.n_basis, false),
    pmat(par.memptr(), bspl.n_basis, bspl.n_basis + 1, false);
  mat pinv_null0 = pinv( bspl.null0 );   // compute pseudoinverse
  spmat.col(0) = pinv_null0 * pmat.col(0);
  spmat.tail_cols(bspl.n_basis - 1) = pinv_null0 * 
    pmat.tail_cols(bspl.n_basis) * pinv_null0.t();
}

/* cond_log_dns: conditional log density
 *
 * return: log density values (double, dim = y.n_elem x x.n_elem) */

mat Item::cond_log_dns(
  vec y,  // y values (double, dim = y.n_elem)
  vec x  // x values (double, dim = x.n_elem)
  )
{
  // initialize
  vec gr;
  mat f(y.n_elem, x.n_elem), he;
  // loop over x and y values
  for (uword j = 0; j < x.n_elem; ++j)
  {
    double log_nc = log_normalize(gr, he, x(j), false);
    for (uword i = 0; i < y.n_elem; ++i)
      f(i, j) = basis_exp(gr, y(i), x(j), false) - log_nc;
  }
  return f;
}
