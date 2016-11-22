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
#include <cmath>
#include <iomanip>


double HermiteInterpolator::getC2C3 (size_t i, double& c2, double& c3) const
{
  if (i+1 > grid.size() || i+1 > values.size() || i+1 > derivs.size())
    return 0.0;

  double h  = grid[i+1] - grid[i];
  double df = (values[i+1] - values[i]) / h;

  c2 = -(2.0*derivs[i] - 3.0*df + derivs[i+1]) / h;
  c3 =  (    derivs[i] - 2.0*df + derivs[i+1]) / (h*h);

  return h;
}


double HermiteInterpolator::evaluate (double x, int derOrder) const
{
  size_t i = 0;
  while (i+2 < grid.size() && x > grid[i+1])
    ++i;

  double c2, c3;
  if (!this->getC2C3(i,c2,c3))
    return 0.0;

  x -= grid[i];
  switch (derOrder) {
  case 0: return values[i] + x*(derivs[i] + x*(c2 + x*c3));
  case 1: return derivs[i] + x*(           2.0*c2 + x*3.0*c3);
  case 2: return                           2.0*c2 + x*6.0*c3;
  case 3: return                                      6.0*c3;
  }

  return 0.0;
}


void HermiteInterpolator::dump (std::ostream& os, size_t n) const
{
  if (n <= grid.size())
    for (size_t i = 0; i < grid.size(); i++)
    {
      std::streamsize oldPrec = os.precision(3);
      std::ios::fmtflags oldF = os.flags(std::ios::right);
      os << grid[i] <<'\t';
      os.precision(15);
      os.flags(std::ios::scientific | std::ios::right);
      os << std::setw(22) << values[i] << std::setw(23) << derivs[i] <<'\n';
      os.precision(oldPrec);
      os.flags(oldF);
    }

  else
  {
    double x, dx = (grid.back() - grid.front())/(n-1);
    for (size_t i = 0; i < n; i++)
    {
      std::streamsize oldPrec = os.precision(3);
      std::ios::fmtflags oldF = os.flags(std::ios::right);
      os << (x = grid.front() + dx*i) <<'\t';
      os.precision(15);
      os.flags(std::ios::scientific | std::ios::right);
      os << std::setw(22) << this->evaluate(x)
         << std::setw(23) << this->evaluateDeriv(x) <<'\n';
      os.precision(oldPrec);
      os.flags(oldF);
    }
  }

  os << std::flush;
}


bool HermiteInterpolator::findMinimum (double& xmin) const
{
#ifdef SP_DEBUG
  std::cout <<"\nHermiteInterpolator::findMinimum:";
#endif

  // Loop to find zeros
  double c2, c3, v, vmin = 1.0e99;
  for (size_t pidx = 0; pidx+1 < grid.size(); pidx++)
  {
    if (this->getC2C3(pidx,c2,c3) == 0.0)
      break;

    // Find (real) roots of quadratic polynomial
    double x = grid[pidx];
    double a = 3.0*c3;
    double b = 2.0*c2 - 6.0*c3*x;
    double c = derivs[pidx] - 2.0*c2*x + 3.0*c3*x*x;
    double d = b*b - 4.0*a*c;
    if (d < 0.0) continue; // complex roots, ignore

    for (int isign = -1; isign <= 1; isign += 2)
    {
      if (a == 0.0 && b != 0.0)
        x = -c/b; // degenerated, b*x + c = 0
      else if (c == 0.0 && a != 0.0)
        x = -b/a; // degenerated, a*x^2 + b*x = 0
      else if (a != 0.0)
        x = 0.5*(isign*sqrt(d) - b)/a;
      else
        break;

#ifdef SP_DEBUG
      std::cout <<" f("<< x <<") = "<< this->evaluate(x);
#endif
      if (x >= grid.front() && x <= grid.back())
        if ((v = this->evaluate(x)) < vmin && this->evaluateDeriv2(x) > 0.0)
        {
          vmin = v;
          xmin = x;
        }
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
