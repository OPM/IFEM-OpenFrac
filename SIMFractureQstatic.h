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
    maxCycle = 50;
    cycleTol = 1.0e-4;
  }

  //! \brief Empty destructor.
  virtual ~SIMFractureQstatic() {}

  //! \brief Parses staggering parameters from an XML element.
  virtual void parseStaggering(const TiXmlElement* elem)
  {
    utl::getAttribute(elem,"tol",cycleTol);
    utl::getAttribute(elem,"max",maxCycle);
  }

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

        return this->checkConvergence(tp,SIM::OK,SIM::CONVERGED) >= SIM::OK;
      }
      else if (this->S1.haveCrackPressure() && rHistory.empty())
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
    else if (status1 != SIM::CONVERGED || status2 != SIM::CONVERGED)
      return SIM::OK;

    int cycle = rHistory.size();

    // Compute residual of the elasticity equation
    this->S1.setMode(SIM::RHS_ONLY);
    if (!this->S1.assembleSystem(tp.time,this->S1.getSolutions(),false))
      return SIM::FAILURE;

    Vector residual;
    if (!this->S1.extractLoadVec(residual))
      return SIM::FAILURE;

    double rNorm1 = residual.norm2();

    // Compute residual of the phase-field equation
    if (!this->S2.setMode(SIM::INT_FORCES))
      return SIM::FAILURE;

    Vectors sol2(1,this->S2.getSolution());
    if (!this->S2.assembleSystem(tp.time,sol2,false))
      return SIM::FAILURE;

    if (!this->S2.extractLoadVec(residual))
      return SIM::FAILURE;

    double rNorm2 = residual.norm2();

    rHistory.push_back(rNorm1+rNorm2);
    double rConv = rHistory.back()/rHistory.front();
    IFEM::cout <<"  cycle="<< cycle <<"  res1="<< rNorm1 <<"  res2="<< rNorm2
               <<"  conv="<< rConv;
    if (cycle > 0)
    {
      double r0 = rHistory.front();
      double r1 = rHistory.back();
      double r2 = rHistory[cycle-1];
      IFEM::cout <<"  beta="<< atan2(cycle*(r2-r1),r0-r1) * 180.0/M_PI;
    }
    IFEM::cout << std::endl;

    if (rConv < cycleTol)
    {
      rHistory.clear();
      return SIM::CONVERGED;
    }
    else if (cycle < maxCycle)
      return SIM::OK;

    std::cerr <<"SIMFractureQstatic::checkConvergence: Did not converge in "
              << maxCycle <<" staggering cycles, bailing.."<< std::endl;
    return SIM::DIVERGED;
  }

private:
  int&      maxCycle; //!< Maximum number of staggering cycles
  double    cycleTol; //!< Residual norm tolerance for the staggering cycles
  RealArray rHistory; //!< Residual norm history for the staggering cycles
};

#endif
