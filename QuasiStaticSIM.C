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
#include "HermiteInterpolator.h"
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
  Vectors tmpSol(1,solution.front()), gNorm;

  if (!model.setMode(SIM::RHS_ONLY))
    return false;

  if (!model.assembleSystem(param.time,tmpSol,false))
    return false;

  if (!model.extractLoadVec(residual))
    return false;

  double fprime0 = residual.dot(linsol);

  if (!model.setMode(SIM::RECOVERY))
    return false;

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

  const size_t numPt = fprime0 > 0.0 ? 21 : 11;
  const double start = fprime0 > 0.0 ? -1.0 : 0.0;
  const double delta = (1.0-start)/(numPt-1);

  RealArray params(numPt), values(numPt), derivs(numPt);
  IFEM::cout <<"\tDoing line search..."<< std::flush;

  for (size_t i = 0; i < numPt; i++)
  {
    params[i] = start + i*delta;
    sol.add(linsol, i == 0 ? start-1.0 : delta);

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
  std::cout <<"\nParameters:\n";
  std::copy(params.begin(),params.end(),
            std::ostream_iterator<double>(std::cout," "));
  std::cout <<"\nValues:\n";
  std::copy(values.begin(),values.end(),
            std::ostream_iterator<double>(std::cout," "));
  std::cout <<"\nDerivatives:\n";
  std::copy(derivs.begin(),derivs.end(),
            std::ostream_iterator<double>(std::cout," "));
  std::cout << std::endl;
#endif
  HermiteInterpolator h(params,values,derivs);
  return h.findMinimum(alpha);
}
