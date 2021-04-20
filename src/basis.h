/* Headers: base and derived class for basis
 *
 * Author: Yang Liu
 *
 * Last modified: 04/14/2021 */

#ifndef BASIS_H
#define BASIS_H

#include "util.h"

// abstract base class Basis
class Basis
{
  public:
    // variables
    uword n_basis;  // number of basis functions
    double lwr, upr;  // lower and upper bounds

    // methods
    virtual rowvec eval(double x) = 0;  // evaluate basis
    virtual Basis *clone() const = 0;  // deep copy

    // constructor and destructor
    Basis() {};
    Basis(uword n_basis_) : n_basis(n_basis_) {};
    Basis(uword n_basis_, double lwr_, double upr_) :
      n_basis(n_basis_), lwr(lwr_), upr(upr_) {};
    virtual ~Basis() {};

};

// derived class Bspline
class Bspline : public Basis
{
  private:
    // variables
    uword order;  // order of spline
    vec knots;  // knot sequence

    // methods
    vec eq_spc_knots();  // equally-spaced knots
    double eval(double x, uword which, 
      uword order);  // single datum, single basis, recursive
    Bspline *clone() const {return new Bspline(*this);}  // deep copy

  public:
    // methods
    rowvec eval(double x);  // evaluate basis
    inline vec get_knots() {return knots;}  // retrieve knots

    // constructor and destructor
    Bspline() {};  // default
    Bspline(uword n_basis_, uword order_, double lwr_, double upr_);
    ~Bspline() {};
};

// derived class Iden
class Iden : public Basis
{
  public:
    // methods
    rowvec eval(double x);  // evaluate basis
    Iden *clone() const {return new Iden(*this);}  // deep copy

    // constructor and destructor
    Iden() {};  // default
    Iden(uword n_basis_);
    ~Iden() {};
};

#endif
