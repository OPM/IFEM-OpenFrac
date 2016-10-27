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
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#if SP_DEBUG > 1
#include <iterator>
#endif


QuasiStaticSIM::QuasiStaticSIM (SIMbase& sim) : NonLinSIM(sim)
{
  numPt = numPtPos = nDump = 0;
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
      utl::getAttribute(child,"nDump",nDump);
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
    if (nDump > 0)
      IFEM::cout <<"\n\tDumping f(alpha) at "<< nDump <<" sampling points.";
  }
  IFEM::cout << std::endl;
}


bool QuasiStaticSIM::evalEnergyFunctional (const TimeDomain& time,
                                           const Vectors& pSol,
                                           double* fVal, double* fDer)
{
  if (fDer)
  {
    if (!model.setMode(SIM::INT_FORCES))
      return false;

    if (!model.assembleSystem(time,pSol,false))
      return false;

    if (!model.extractLoadVec(residual))
      return false;

    *fDer = -residual.dot(linsol);
  }

  if (fVal)
  {
    if (!model.setMode(SIM::RECOVERY))
      return false;

    Vectors gNorm;
    if (!model.solutionNorms(time,pSol,gNorm))
      return false;

    *fVal = gNorm.front()(1) + gNorm.front()(6);
  }

  return true;
}


bool QuasiStaticSIM::lineSearch (TimeStep& param)
{
  alpha = 1.0;
  if (numPtPos < 2)
    return true; // No line search

  // Make a temporary copy of the primary solution
  Vectors tmpSol(1,solution.front());
  Vector& solVec = tmpSol.front();
  double curr, prev;

  // Evaluate f(alpha) at alpha=1.0
  solVec.add(linsol);
  if (!this->evalEnergyFunctional(param.time,tmpSol,&curr))
    return false;

  // Evaluate f(alpha) at alpha=0.0
  solVec.add(linsol,-1.0);
  if (!this->evalEnergyFunctional(param.time,tmpSol,&prev))
    return false;

#ifdef SP_DEBUG
  std::cout <<"\tLine search? curr: "<< curr <<" prev: "<< prev << std::endl;
#endif
  if (curr < prev)
    return true; // No line search needed in this iteration

  // Evaluate f'(alpha) at alpha=0.0
  if (!this->evalEnergyFunctional(param.time,tmpSol,nullptr,&curr))
    return false;

  IFEM::cout <<"\tDoing line search, f'(0) = "<< curr <<" :"<< std::flush;

  // If f'(alpha) < 0 then use sampling interval [0,1] else use [-1,1]
  const size_t nVal = curr < 0.0 ? numPtPos : numPt;
  RealArray prm(params.begin()+numPt-nVal,params.end());
  RealArray values(nVal), derivs(nVal);

  // Evaluate f(alpha) and f'(alpha) forall alpha in [prm.front(),prm.back()]
  for (size_t i = 0; i < nVal; i++)
  {
    if (i > 0)
      solVec.add(linsol,prm[i] - prm[i-1]);
    else if (prm.front() != 0.0)
      solVec.add(linsol,prm.front());
    if (!this->evalEnergyFunctional(param.time,tmpSol,&values[i],&derivs[i]))
      return false;
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
  if (nDump > 0)
  {
    // Dump the scalar energy functional to file for plotting
    char filename[64];
    sprintf(filename,"f_alpha_s%d_i%d.dat",param.step,param.iter);
    std::ofstream os(filename);
    h.dump(os,nDump);
  }
  return h.findMinimum(alpha);
}
