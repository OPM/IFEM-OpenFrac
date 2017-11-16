// $Id$
//==============================================================================
//!
//! \file SIMFractureQstatic.h
//!
//! \date Aug 9 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for staggered quasti-static brittle fracture problems.
//!
//==============================================================================

#ifndef _SIM_FRACTURE_QSTATIC_H_
#define _SIM_FRACTURE_QSTATIC_H_

#include "SIMFractureDynamics.h"
#include "SIMCoupledSI.h"


/*!
  \brief Driver class for staggered quasi-static fracture simulators.
*/

template<class SolidSlv, class PhaseSlv>
class SIMFractureQstatic : public SIMFracture<SolidSlv,PhaseSlv,SIMCoupledSI>
{
  //! Convenience type
  typedef SIMFracture<SolidSlv,PhaseSlv,SIMCoupledSI> CoupledSIM;

public:
  //! \brief The constructor forwards to the parent class contructor.
  SIMFractureQstatic(SolidSlv& s1, PhaseSlv& s2, const std::string& input)
    : CoupledSIM(s1,s2,input), maxCycle(this->maxIter)
  {
    maxCycle = maxIt = 50;
    cycleTol = 1.0e-4;
  }

  //! \brief Empty destructor.
  virtual ~SIMFractureQstatic() {}

  //! \brief Parses staggering parameters from an XML element.
  virtual void parseStaggering(const TiXmlElement* elem)
  {
    this->CoupledSIM::parseStaggering(elem);
    utl::getAttribute(elem,"tol",cycleTol);
    utl::getAttribute(elem,"max",maxCycle);
    maxIt = maxCycle;
  }

  //! \brief Enable/disable the staggering iteration cycles.
  virtual void enableStaggering(bool enable) { maxCycle = enable ? maxIt : 0; }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool firstS1 = true)
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

    return this->SIMCoupledSI<SolidSlv,PhaseSlv>::solveStep(tp,firstS1);
  }

  //! \brief Checks if the coupled simulator has converged.
  virtual SIM::ConvStatus checkConvergence(const TimeStep& tp,
                                           SIM::ConvStatus status1,
                                           SIM::ConvStatus status2)
  {
    if (status1 == SIM::FAILURE || status2 == SIM::FAILURE)
      return SIM::FAILURE;
    else if (status1 == SIM::DIVERGED || status2 == SIM::DIVERGED)
      return SIM::DIVERGED;

    //! \brief Calculate and print solution and residual norms
    double conv = this->calcResidual(tp,true);
    if (conv < 0.0)
      return SIM::FAILURE;
    else if (conv < fabs(cycleTol))
      return SIM::CONVERGED;
    else if (tp.iter < maxCycle)
      return SIM::OK;
    else if (cycleTol < 0.0 || maxCycle == 0)
      return SIM::CONVERGED; // Continue after maximum number of cycles

    std::cerr <<"SIMFractureQstatic::checkConvergence: Did not converge in "
              << maxCycle <<" staggering cycles, bailing.."<< std::endl;
    return SIM::DIVERGED;
  }

private:
  int&   maxCycle; //!< Maximum number of staggering cycles
  int    maxIt;    //!< Copy of \a maxCycle
  double cycleTol; //!< Residual norm tolerance for the staggering cycles
  double E0;       //!< Energy norm of initial staggering cycle
  double Ec;       //!< Energy norm of current staggering cycle
  double Ep;       //!< Energy norm of previous staggering cycle
};


/*!
  \brief Driver class for quasi-static fracture simulators.
*/

template<class SolidSlv, class PhaseSlv>
class SIMFractureMiehe : public SIMFracture<SolidSlv,PhaseSlv,SIMCoupled>
{
  //! Convenience type
  typedef SIMFracture<SolidSlv,PhaseSlv,SIMCoupled> CoupledSIM;

public:
  //! \brief The constructor forwards to the parent class contructor.
  SIMFractureMiehe(SolidSlv& s1, PhaseSlv& s2, const std::string& input)
    : CoupledSIM(s1,s2,input)
  {
    doStg = true;
    numCycle = 2;
    cycleTol = 1.0e-4;
  }

  //! \brief Empty destructor.
  virtual ~SIMFractureMiehe() {}

  //! \brief Parses staggering parameters from an XML element.
  virtual void parseStaggering(const TiXmlElement* elem)
  {
    this->CoupledSIM::parseStaggering(elem);
    utl::getAttribute(elem,"tol",cycleTol);
    utl::getAttribute(elem,"max",numCycle);
    if (cycleTol < 0.0) cycleTol = -cycleTol;
  }

  //! \brief Enable/disable the staggering iteration cycles.
  virtual void enableStaggering(bool enable = true) { doStg = enable; }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool = true)
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

    return cycleTol == 0.0 ? this->calcResidual(tp) >= 0.0 : true;
  }

private:
  bool   doStg;    //!< Toggle for enabling/disabling staggering cycles
  int    numCycle; //!< Number of staggering cycles
  double cycleTol; //!< Residual norm tolerance for the staggering cycles
};

#endif
