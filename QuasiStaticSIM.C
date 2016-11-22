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
#include "FractureElasticityVoigt.h"
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


QuasiStaticSIM::QuasiStaticSIM (SIMbase& sim) : NonLinSIM(sim,L2)
{
  numPt = numPtPos = nDump = 0;
  version = 1;
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
      utl::getAttribute(child,"version",version);
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
    if (version == 1)
      IFEM::cout <<"\n\tEvaluating f'(alpha) as {F}_int * {u}_corr.";
    else if (version == 2)
      IFEM::cout <<"\n\tEvaluating f'(alpha) by direct integration.";
  }
  IFEM::cout << std::endl;
}


bool QuasiStaticSIM::evalEnergyFunctional (const TimeDomain& time,
                                           const Vectors& pSol,
                                           double* fVal, double* fDer)
{
  if (fDer && version == 1)
  {
    // Integrate f'(alpha) as the dot-product between internal force vector
    // and the solution correction vector
    if (!model.setMode(SIM::INT_FORCES))
      return false;

    if (!model.assembleSystem(time,pSol,false))
      return false;

    Vector forces;
    if (!model.extractLoadVec(forces))
      return false;

    *fDer = -forces.dot(linsol);
  }

  if (fVal || version == 2)
  {
    FractureElasticNorm::dirDer = fDer && version == 2;

    if (!model.setMode(SIM::RECOVERY))
      return false;

    Vectors gNorm;
    if (!model.solutionNorms(time,pSol,gNorm))
      return false;

    if (fVal)
      *fVal = gNorm.front()(1) + gNorm.front()(6);
    if (fDer && version == 2)
      *fDer = gNorm.front()(5); // f'(alpha) is integrated by the norm class

    FractureElasticNorm::dirDer = false;
  }

  return true;
}


bool QuasiStaticSIM::lineSearch (TimeStep& param)
{
  alpha = alphaO = 1.0;
  if (numPtPos < 2)
    return true; // No line search

  FractureElasticNorm::extEnr = false; // Disable external energy calculation,
  // since it uses a path integral valid at the converged configuration only

  // Make a temporary copy of the primary solution
  Vectors tmpSol(version == 2 ? 2 : 1, solution.front());
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
  {
    FractureElasticNorm::extEnr = true; // Enable external energy calculation
    return true; // No line search needed in this iteration
  }

  if (version == 2) // Store the solution correction in tmpSol
    tmpSol.back() = linsol;

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
  FractureElasticNorm::extEnr = true; // Enable external energy calculation

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

  // Find the value of alpha that minimizes the function f(alpha)
  if (h.findMinimum(alpha))
    alphaO = alpha;
  else
    return false;

  // Evaluate the residual at the new solution for the convergence check
  if (!model.setMode(SIM::RHS_ONLY))
    return false;

  solVec.add(linsol,alpha-prm.back());
  if (!model.assembleSystem(param.time,tmpSol,false))
    return false;

  return model.extractLoadVec(residual);
}
