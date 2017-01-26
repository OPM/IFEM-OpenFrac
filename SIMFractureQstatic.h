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
    E0 = Ec = Ep = 0.0;
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

    // Compute residual of the elasticity equation
    this->S1.setMode(SIM::RHS_ONLY);
    if (!this->S1.assembleSystem(tp.time,this->S1.getSolutions(),false))
      return SIM::FAILURE;

    Vector residual;
    if (!this->S1.extractLoadVec(residual))
      return SIM::FAILURE;

    double rNorm1 = residual.norm2();
    double eNorm1 = this->S1.extractScalar();

    // Compute residual of the phase-field equation
    if (!this->S2.setMode(SIM::INT_FORCES))
      return SIM::FAILURE;

    Vectors sol2(1,this->S2.getSolution());
    if (!this->S2.assembleSystem(tp.time,sol2,false))
      return SIM::FAILURE;

    if (!this->S2.extractLoadVec(residual))
      return SIM::FAILURE;

    double rNorm2 = residual.norm2();
    double eNorm2 = this->S2.extractScalar();

    double rConv = rNorm1 + rNorm2;
    double eConv = eNorm1 + eNorm2;
    IFEM::cout <<"  cycle "<< tp.iter <<": Res = "<< rNorm1 <<" + "<< rNorm2
               <<" = "<< rConv <<"  E = "<< eNorm1 <<" + "<< eNorm2
               <<" = "<< eConv;
    if (tp.iter == 0)
      E0 = eConv;
    else
    {
      Ep = tp.iter > 1 ? Ec : E0;
      Ec = eConv;
      IFEM::cout <<"  beta="<< atan2(tp.iter*(Ep-Ec),E0-Ec) * 180.0/M_PI;
    }
    IFEM::cout << std::endl;

    if (rConv < cycleTol)
      return SIM::CONVERGED;
    else if (tp.iter < maxCycle)
      return SIM::OK;

    std::cerr <<"SIMFractureQstatic::checkConvergence: Did not converge in "
              << maxCycle <<" staggering cycles, bailing.."<< std::endl;
    return SIM::DIVERGED;
  }

private:
  int&   maxCycle; //!< Maximum number of staggering cycles
  double cycleTol; //!< Residual norm tolerance for the staggering cycles
  double E0;       //!< Energy norm of initial staggering cycle
  double Ec;       //!< Energy norm of current staggering cycle
  double Ep;       //!< Energy norm of previous staggering cycle
};

#endif
