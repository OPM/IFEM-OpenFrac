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
  SIMFractureQstatic(SolidSlv& s1, PhaseSlv& s2, const std::string& input);

  //! \brief Empty destructor.
  virtual ~SIMFractureQstatic() {}

  //! \brief Parses staggering parameters from an XML element.
  virtual void parseStaggering(const TiXmlElement* elem);

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool firstS1 = true);

  //! \brief Checks if the coupled simulator has converged.
  virtual SIM::ConvStatus checkConvergence(const TimeStep& tp,
                                           SIM::ConvStatus status1,
                                           SIM::ConvStatus status2);

private:
  int    maxInc;   //!< Max number of cycles with increasing convergence norm
  double lastConv; //!< Previous convergence norm value
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
  SIMFractureMiehe(SolidSlv& s1, PhaseSlv& s2, const std::string& input);

  //! \brief Empty destructor.
  virtual ~SIMFractureMiehe() {}

  //! \brief Parses staggering parameters from an XML element.
  virtual void parseStaggering(const TiXmlElement* elem);

  //! \brief Enable/disable the staggering iteration cycles.
  virtual void enableStaggering(bool enable = true) { doStg = enable; }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool = true);

private:
  bool   doStg;    //!< Toggle for enabling/disabling staggering cycles
  int    numCycle; //!< Number of staggering cycles
  double cycleTol; //!< Residual norm tolerance for the staggering cycles
};

#endif
