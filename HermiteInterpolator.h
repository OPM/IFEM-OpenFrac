// $Id$
//==============================================================================
//!
//! \file HermiteInterpolator.h
//!
//! \date Oct 18 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Piece-wise cubic hermite interpolation utilities.
//!
//==============================================================================

#ifndef HERMITE_INTERPOLATOR_H_
#define HERMITE_INTERPOLATOR_H_

#include <vector>


/*!
  \brief Piece-wise cubic hermite interpolation class.
*/

class HermiteInterpolator
{
public:
  //! \brief The constructor initializes the data members.
  HermiteInterpolator(const std::vector<double>& grid_,
                      const std::vector<double>& values_,
                      const std::vector<double>& derivs_)
    : grid(grid_), values(values_), derivs(derivs_) {}

  //! \brief Evaluates interpolation polynomial at \a x.
  double evaluate(double x) const { return this->evaluate(x,0); }
  //! \brief Evaluates derivative of interpolation polynomial at \a x.
  double evaluateDeriv(double x) const { return this->evaluate(x,1); }
  //! \brief Evaluates 2nd derivative of interpolation polynomial at \a x.
  double evaluateDeriv2(double x) const { return this->evaluate(x,2); }

  //! \brief Finds the minimum of a cubic hermite interpolant.
  bool findMinimum(double& xmin) const;

private:
  //! \brief Evaluates derOrder'th derivative of the interpolation polynomial.
  double evaluate(double x, int derOrder) const;

protected:
  std::vector<double> grid;   //!< Grid points
  std::vector<double> values; //!< Function values
  std::vector<double> derivs; //!< Function derivatives
};

#endif
