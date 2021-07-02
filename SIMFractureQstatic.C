// $Id$
//==============================================================================
//!
//! \file SIMFractureQstatic.C
//!
//! \date Aug 9 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for staggered quasti-static brittle fracture problems.
//!
//==============================================================================

#include "SIMFractureQstatic.h"

#include "SIMDynElasticity.h"
#include "SIMPhaseField.h"
#include "SIMExplPhaseField.h"
#include "SIMPoroElasticity.h"

#include "GenAlphaSIM.h"
#include "HHTSIM.h"
#include "NewmarkNLSIM.h"
#include "NonLinSIM.h"
#include "SIM2D.h"
#include "SIM3D.h"


template<class SolidSlv, class PhaseSlv>
SIMFractureQstatic<SolidSlv,PhaseSlv>::
SIMFractureQstatic (SolidSlv& s1, PhaseSlv& s2, const std::string& input)
  : CoupledSIM(s1,s2,input)
{
  maxInc   = 0;
  cycleTol = 1.0e-4;
}


template<class SolidSlv, class PhaseSlv>
void SIMFractureQstatic<SolidSlv,PhaseSlv>::
parseStaggering (const TiXmlElement* elem)
{
  this->CoupledSIM::parseStaggering(elem);
  this->CoupledSIM::parseIterations(elem);
  utl::getAttribute(elem,"maxInc",maxInc);
  utl::getAttribute(elem,"tol",cycleTol);
}


template<class SolidSlv, class PhaseSlv>
bool SIMFractureQstatic<SolidSlv,PhaseSlv>::
solveStep (TimeStep& tp, bool firstS1)
{
  if (tp.step == 1)
  {
    // Only solve the elasticity problem in the first step,
    if (this->S2.hasIC("phasefield")) // if an initial phase field is given
    {
      IFEM::cout <<"\n  Initial phase field..."<< std::endl;
      if (!this->S2.postSolve(tp))
        return false;

      TimeStep myTp(tp); // Make a copy to avoid changing the cycle counter
      if (!this->S1.solveStep(myTp))
        return false;

      if (this->calcResidual(tp) < 0.0)
        return false;

      tp.time.first = false;
      return true;
    }
    else if (this->S1.haveCrackPressure())
      // Start the initial step by solving the phase-field first
      if (!this->S2.solveStep(tp,false))
        return false;
  }
  else // solve the phase-field equation first, if an initial field is given
    firstS1 = !this->S2.hasIC("phasefield");

  if (!this->SIMCoupledSI<SolidSlv,PhaseSlv>::solveStep(tp,firstS1))
    return false;

  CoupledSIM::doStop = this->S2.checkStopCriterion();
  return true;
}


template<class SolidSlv, class PhaseSlv>
SIM::ConvStatus SIMFractureQstatic<SolidSlv,PhaseSlv>::
checkConvergence (const TimeStep& tp,
                  SIM::ConvStatus status1,
                  SIM::ConvStatus status2)
{
  if (status1 == SIM::FAILURE || status2 == SIM::FAILURE)
    return SIM::FAILURE;
  else if (status1 == SIM::DIVERGED || status2 == SIM::DIVERGED)
    return SIM::DIVERGED;

  double conv = this->calcResidual(tp,true);
  if (conv < 0.0)
    return SIM::FAILURE;
  else if (conv < fabs(cycleTol))
    return SIM::CONVERGED;

  static int numIncr = 0;
  if (tp.iter == 0 || conv <= lastConv)
    numIncr = 0;
  else
    numIncr++;

  lastConv = conv;
  if (numIncr > maxInc && maxInc > 0)
  {
    IFEM::cout <<"  ** The residual increases in more than "<< maxInc
               <<" cycles, giving up and continuing..."<< std::endl;
    return SIM::CONVERGED; // Abort cycles and continue with next step
  }

  int maxCycle = this->getMaxit(tp.step);
  if (tp.iter < maxCycle)
    return SIM::OK; // Continue with next cycle
  else if (cycleTol < 0.0 || maxCycle == 0 || maxInc > 0)
    return SIM::CONVERGED; // Continue after maximum number of cycles

  std::cerr <<"SIMFractureQstatic::checkConvergence: Did not converge in "
            << maxCycle <<" staggering cycles, bailing.."<< std::endl;
  return SIM::DIVERGED;
}


template<class SolidSlv, class PhaseSlv>
SIMFractureMiehe<SolidSlv,PhaseSlv>::
SIMFractureMiehe (SolidSlv& s1, PhaseSlv& s2, const std::string& input)
  : CoupledSIM(s1,s2,input)
{
  doStg = true;
  numCycle = 2;
  cycleTol = 1.0e-4;
}


template<class SolidSlv, class PhaseSlv>
void SIMFractureMiehe<SolidSlv,PhaseSlv>::
parseStaggering (const TiXmlElement* elem)
{
  this->CoupledSIM::parseStaggering(elem);
  utl::getAttribute(elem,"tol",cycleTol);
  utl::getAttribute(elem,"max",numCycle);
  if (cycleTol < 0.0) cycleTol = -cycleTol;
}


template<class SolidSlv, class PhaseSlv>
bool SIMFractureMiehe<SolidSlv,PhaseSlv>::
solveStep (TimeStep& tp, bool)
{
  if (tp.step == 1)
  {
    // Only solve the elasticity problem in the first step,
    if (this->S2.hasIC("phasefield")) // if an initial phase field is given
    {
      IFEM::cout <<"\n  Initial phase field..."<< std::endl;
      if (!this->S2.postSolve(tp))
        return false;

      TimeStep myTp(tp); // Make a copy to avoid changing the cycle counter
      if (!this->S1.solveStep(myTp))
        return false;

      if (cycleTol > 0.0 && this->calcResidual(tp) < 0.0)
        return false;

      tp.time.first = false;
      return true;
    }
    else if (this->S1.haveCrackPressure())
      // Start the initial step by solving the phase-field first
      if (!this->S2.solveStep(tp,false))
        return false;
  }

  tp.iter = 0; // Solve the predictor step for the elasticity problem
  if (this->S1.solveIteration(tp,1) <= SIM::DIVERGED)
    return false;

  // Update strain energy density for the predictor step
  if (!this->S1.updateStrainEnergyDensity(tp))
    return false;

  // Solve the phase-field problem
  if (!this->S2.solveStep(tp,false))
    return false;

  TimeStep myTp(tp); // Make a copy to avoid changing the cycle counter
  ++myTp.iter; // Iterate the elasticity problem (corrector steps)
  if (this->S1.solveIteration(myTp,2) <= SIM::DIVERGED)
    return false;

  double conv = cycleTol > 0.0 ? this->calcResidual(tp,true) : 1.0;
  for (tp.iter = 1; tp.iter < numCycle && conv > cycleTol && doStg; tp.iter++)
  {
    // Solve the phase-field problem
    if (!this->S2.solveStep(tp,false))
      return false;

    // Solve the elasticity problem
    if (this->S1.solveIteration(tp,3) <= SIM::DIVERGED)
      return false;

    // Check if the staggering cycles have converged
    if (cycleTol > 0.0)
      conv = this->calcResidual(tp,true);
  }
  if (conv < 0.0) return false;

  tp.time.first = false;
  this->S1.postSolve(tp);
  this->S2.postSolve(tp);
  if (cycleTol == 0.0 && this->calcResidual(tp) < 0.0)
    return false;

  CoupledSIM::doStop = this->S2.checkStopCriterion();
  return true;
}


//! \brief Helper macro doing the actual instantation.
#define INSTANCE_FULL(TYPE,DIM,SIM,ELSIM) \
  template class TYPE<SIMDynElasticity<DIM,SIM,ELSIM<DIM>>, \
                      SIMPhaseField<DIM>>; \
  template class TYPE<SIMDynElasticity<DIM,SIM,ELSIM<DIM>>, \
                      SIMExplPhaseField>;

//! \brief Helper macro adding dimensionality.
#define INSTANCE_DIM(TYPE,SIM,ELSIM) \
  INSTANCE_FULL(TYPE,SIM2D,SIM,ELSIM) \
  INSTANCE_FULL(TYPE,SIM3D,SIM,ELSIM)

//! \brief Helper macro adding nonlinear solver type.
#define INSTANCE_TYPE(TYPE,ELSIM) \
  INSTANCE_DIM(TYPE,GenAlphaSIM,ELSIM) \
  INSTANCE_DIM(TYPE,HHTSIM,ELSIM) \
  INSTANCE_DIM(TYPE,NewmarkSIM,ELSIM) \
  INSTANCE_DIM(TYPE,NewmarkNLSIM,ELSIM) \
  INSTANCE_DIM(TYPE,NonLinSIM,ELSIM)

//! \brief Helper macro to instance for a given elasticity type.
#define INSTANCE(ELSIM) \
  INSTANCE_TYPE(SIMFractureQstatic,ELSIM) \
  INSTANCE_TYPE(SIMFractureMiehe,ELSIM)

INSTANCE(SIMElasticityWrap)
INSTANCE(SIMPoroElasticity)
