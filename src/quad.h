/* Headers: class _Quad_
 *
 * Author: Yang Liu
 *
 * Last modified: 08/26/2020 */

#ifndef QUAD_H
#define QUAD_H

#include <RcppArmadillo.h>
using namespace std;
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// class _Quad_
class Quad
{
  private:
    // variables
    uword n_quad;  // number of quadrature points

    // functions
    void gauss_legendre(double lower, double upper);  // Gauss-Legendre quadrature

  public:
    // constructors
    Quad(){};  // default
    Quad(uword n_quad_, double lower, double upper);

    // variables
    vec node;   // quadrature nodes
    vec weight;  // quadrature weights

  friend class Item;
  friend class Test;
};

#endif
