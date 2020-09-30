/* Interfaces for internal classes
 *
 * Author: Yang Liu
 *
 * Last modified: 08/29/2020 */

#include "test.h"

/* score: compute EAP score
 *
 * update: x (double, dim = n_obsn) */

vec Test::score()
{
  vec x = estep_wt.t() * quad.node;
  return x;
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

/* cond_dns: conditional density
 *
 * return: conditional density values */

mat Item::cond_dns(
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
    {
      f(i, j) = std::exp(basis_exp(gr, y(i), x(j), false) - 
        log_nc);
    }
  }
  return f;
}
