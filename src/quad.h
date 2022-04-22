/* Headers: class Quad
 *
 * Author: Yang Liu
 *
 * Last modified: 04/21/2022 */

#ifndef QUAD_H
#define QUAD_H

#include "util.h"

// abstract base class Quad
class Quad
{
  public:
    // variables
    uword n_quad1, n_dim, 
      n_quad; // # of quad points per dimension, dimension, total # of quad points 
    double lwr, upr;  // upper and lower limit
    mat node;   // quadrature nodes
    vec weight;  // quadrature weights

    // method
    virtual Quad *clone() const = 0;  // deep copy

    // constructor and destructor
    Quad() {};  // default
    Quad(uword n_quad1_, uword n_dim_);
    Quad(uword n_quad1_) : Quad(n_quad1_, 1) {};
    Quad(uword n_quad1_, uword n_dim_, double lwr_, double upr_);
    virtual ~Quad() {};
};

// Gauss-Legendre quadrature
class GaussLegendre : public Quad
{
  public:
    // method
    GaussLegendre *clone() const {return new GaussLegendre(*this);}  // deep copy
    
    // constructor and destructor
    GaussLegendre() {};  // default
    GaussLegendre(uword n_quad1_, uword n_dim_, double lwr_, double upr_);
    virtual ~GaussLegendre() {};
};

// constant quadrature (unnormalized)
class Const : public Quad
{
  public:
    // method
    Const *clone() const {return new Const(*this);}  // deep copy
    
    // constructor and destructor
    Const() {};  // default
    Const(uword n_quad1_);
    virtual ~Const() {};
};

#endif
