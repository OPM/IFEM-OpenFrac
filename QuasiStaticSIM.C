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
#include "Utilities.h"
#include "tinyxml.h"
#include "IFEM.h"
#include <cstring>
#include <cstdlib>
#if SP_DEBUG > 1
#include <iterator>
#endif


QuasiStaticSIM::QuasiStaticSIM (SIMbase& sim) : NonLinSIM(sim)
{
  numPt = numPtPos = 0;
}


bool QuasiStaticSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),inputContext))
    return model.parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"linesearch"))
    {
      if (utl::getAttribute(child,"numPt",numPt) && numPt > 2)
      {
        // Using a uniform point distribution over the domain [-1,1]
        double delta = 2.0/(numPt-1);
        for (size_t i = 0; i < numPt; i++)
          params.push_back(i*delta - 1.0);
        numPtPos = numPt%2 + numPt/2;
      }
      else if (child->FirstChild())
      {
        // User-defined non-uniform point distribution
        std::string value(child->FirstChild()->Value());
        char* cval = strtok(const_cast<char*>(value.c_str())," ");
        for (int i = 0; cval; i++, cval = strtok(nullptr," "))
          if (params.empty())
            params.push_back(atof(cval));
          else if ((alpha = atof(cval)) > params.back())
          {
            params.push_back(alpha);
            if (alpha >= 0.0) numPtPos++;
          }
        numPt = params.size();
      }
    }

  return this->NonLinSIM::parse(elem);
}


void QuasiStaticSIM::printProblem () const
{
  model.printProblem();

  IFEM::cout <<"Quasi-static nonlinear solution driver";
  if (numPtPos > 1)
  {
    IFEM::cout <<" with line search\n\tParameters: "<< params.front();
    for (size_t i = 1; i < params.size(); i++)
      IFEM::cout << (i%10 ? " " : "\n\t            ") << params[i];
  }
  IFEM::cout << std::endl;
}


bool QuasiStaticSIM::lineSearch (TimeStep& param)
{
  alpha = 1.0;
  if (numPtPos < 2)
    return true; // No line search

  Vectors tmpSol(1,solution.front()), gNorm;

  if (!model.setMode(SIM::RHS_ONLY))
    return false;

  if (!model.assembleSystem(param.time,tmpSol,false))
    return false;

  if (!model.extractLoadVec(residual))
    return false;

  const size_t nVal = residual.dot(linsol) > 0.0 ? numPt : numPtPos;

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
  if (curr < prev)
    return true; // No line search needed in this iteration

  RealArray prm(params.begin()+numPt-nVal,params.end());
  RealArray values(nVal), derivs(nVal);
  IFEM::cout <<"\tDoing line search..."<< std::flush;

  for (size_t i = 0; i < nVal; i++)
  {
    sol.add(linsol, i == 0 ? prm.front()-1.0 : prm[i]-prm[i-1]);

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
  std::copy(prm.begin(),prm.end(),
            std::ostream_iterator<double>(std::cout," "));
  std::cout <<"\nValues:\n";
  std::copy(values.begin(),values.end(),
            std::ostream_iterator<double>(std::cout," "));
  std::cout <<"\nDerivatives:\n";
  std::copy(derivs.begin(),derivs.end(),
            std::ostream_iterator<double>(std::cout," "));
  std::cout << std::endl;
#endif
  HermiteInterpolator h(prm,values,derivs);
  return h.findMinimum(alpha);
}
