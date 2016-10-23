// $Id$
//==============================================================================
//!
//! \file HermiteInterpolator.C
//!
//! \date Oct 18 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Piece-wise cubic hermite interpolation utilities.
//!
//==============================================================================

#include "HermiteInterpolator.h"
#include "IFEM.h"


double HermiteInterpolator::evaluate (double x, int derOrder) const
{
  size_t i = 0;
  while (i+2 < grid.size() && x > grid[i+1])
    ++i;

  double h  = grid[i+1] - grid[i];
  double df = (values[i+1] - values[i]) / h;

  double c2 = -(2.0*derivs[i] - 3.0*df + derivs[i+1]) / h;
  double c3 =  (    derivs[i] - 2.0*df + derivs[i+1]) / (h*h);

  x -= grid[i];
  switch (derOrder) {
  case 0: return values[i] + x*(derivs[i] + x*(c2 + x*c3));
  case 1: return derivs[i] + x*(           2.0*c2 + x*3.0*c3);
  case 2: return                           2.0*c2 + x*6.0*c3;
  case 3: return                                      6.0*c3;
  }

  return 0.0;
}


bool HermiteInterpolator::findMinimum (double& xmin, unsigned int maxIts,
                                       double absTol, double relTol) const
{
  const double epsZero = 1.0e-10;

#ifdef SP_DEBUG
  std::cout <<"\nHermiteInterpolator::findMinimum:";
#endif

  // Newton loop to find zeros
  double v, vmin = 1.0e99;
  for (size_t pidx = 0; pidx < grid.size(); pidx++)
  {
    double x = grid[pidx], dx = 1.0;
    for (unsigned int its = 0; its < maxIts; its++)
      if ((dx > absTol || dx < -absTol) && (dx > relTol*x || dx < -relTol*x))
      {
        double I  = this->evaluateDeriv(x);
        double I2 = this->evaluateDeriv2(x);
        if (I2 < epsZero && I2 > -epsZero)
          return false; // breakdown - probably a linear function
        dx = I / I2;
        x -= dx;
      }
      else // converged
      {
#ifdef SP_DEBUG
        std::cout <<" f("<< x <<") = "<< this->evaluate(x);
#endif
        if (x >= grid.front() && x <= grid.back())
          if ((v = this->evaluate(x)) < vmin && this->evaluateDeriv2(x) > 0.0)
          {
            vmin = v;
            xmin = x;
          }
        break;
      }
  }

  if (vmin < 1.0e99)
  {
#ifndef SP_DEBUG
    IFEM::cout <<"  alpha^* = "<< xmin << std::endl;
#endif
    return true;
  }

  IFEM::cout <<"\n\t*** No energy minimum found in [u";
  if (grid.front() == 0.0)
    IFEM::cout <<",u";
  else if (grid.front() == -1.0)
    IFEM::cout <<"-du,u";
  else
    IFEM::cout << grid.front() <<"*du,u";
  if (grid.back() == 1.0)
    IFEM::cout <<"+du]";
  else
    IFEM::cout << grid.front() <<"*du]";
  IFEM::cout << std::endl;
  return false;
}
