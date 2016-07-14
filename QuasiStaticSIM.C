// $Id$
//==============================================================================
//!
//! \file QuasiStaticSIM.C
//!
//! \date Sep 22 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Quasi-static solution driver for fracture elasticity simulators.
//!
//==============================================================================

#include "QuasiStaticSIM.h"
#include "CubicMinimum.h"
#include "SIMoutput.h"
#include "TimeStep.h"
#include "IFEM.h"
#if SP_DEBUG > 1
#include <iterator>
#endif


QuasiStaticSIM::QuasiStaticSIM (SIMbase& sim) : NonLinSIM(sim)
{
  eta = 1.0;
}


void QuasiStaticSIM::printProblem () const
{
  model.printProblem();

  IFEM::cout <<"Quasi-static nonlinear solution driver";
  if (eta > 0.0) IFEM::cout <<" with line search";
  IFEM::cout << std::endl;
}


bool QuasiStaticSIM::lineSearch (TimeStep& param)
{
  if (eta <= 0.0)
    return true; // No line search

  alpha = 1.0;

  if (!model.setMode(SIM::RECOVERY))
    return false;

  Vectors tmpSol(1,solution.front()), gNorm;
  if (!model.solutionNorms(param.time,tmpSol,gNorm))
    return false;

  Vector& sol  = tmpSol.front();
  Vector& norm = gNorm.front();
  double  prev = norm(1) + norm(6);

  sol.add(linsol);
  if (!model.solutionNorms(param.time,tmpSol,gNorm))
    return false;

  double curr = norm(1) + norm(6);

#ifdef SP_DEBUG
  std::cout <<"\tLine search? curr: "<< curr <<" prev: "<< prev << std::endl;
#endif
  if (param.iter < 2 || curr < prev)
    return true; // No line search needed in this iteration
#ifdef SP_DEBUG
  std::cout <<"\tDoing line search."<< std::endl;
#endif

  const size_t numPt = 10;
  const double delta = 2.0/(numPt-1);

  RealArray params(numPt), values(numPt), derivs(numPt);

  for (size_t i = 0; i < numPt; i++)
  {
    params[i] = i*delta - 1.0;
    sol.add(linsol, i == 0 ? -2.0 : delta);

    if (!model.setMode(SIM::RHS_ONLY))
      return false;

    if (!model.assembleSystem(param.time,tmpSol,false))
      return false;

    if (!model.extractLoadVec(residual))
      return false;

    if (!model.setMode(SIM::RECOVERY))
      return false;

    if (!model.solutionNorms(param.time,Vectors(1,sol),gNorm))
      return false;

    values[i] = norm(1) + norm(6);
    derivs[i] = residual.dot(linsol);
  }

#if SP_DEBUG > 1
  std::cout <<"\nValues:\n";
  std::copy(values.begin(),values.end(),
            std::ostream_iterator<double>(std::cout," "));
  std::cout <<"\nDerivatives:\n";
  std::copy(derivs.begin(),derivs.end(),
            std::ostream_iterator<double>(std::cout," "));
  std::cout << std::endl;
#endif

  return CubicMinimum::Find(alpha,params,values,derivs);
}
