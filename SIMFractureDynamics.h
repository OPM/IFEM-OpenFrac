// $Id$
//==============================================================================
//!
//! \file SIMFactureDynamics.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for fracture-dynamic problems.
//!
//==============================================================================

#ifndef _SIM_FRACTURE_DYNAMICS_H_
#define _SIM_FRACTURE_DYNAMICS_H_

#include "SIMCoupled.h"


/*!
  \brief Driver class for fracture dynamics simulators.
  \details A fracture dynamics simulator is a coupling between a phase field solver
  and an elasticity solver.
*/

template<class SolidSolver, class PhaseSolver>
class SIMFractureDynamics : public SIMCoupled<SolidSolver,PhaseSolver>
{
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMFractureDynamics(SolidSolver& s1, PhaseSolver& s2)
    : SIMCoupled<SolidSolver,PhaseSolver>(s1,s2) {}

  //! \brief Empty destructor.
  virtual ~SIMFractureDynamics() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies()
  {
    this->S2.registerDependency(&this->S1, "tensile", 1);
    // This is defined on integration point and not on control points.
    // It is a global vector across all patches on the process.
    // Use an explicit call instead of normal couplings for this.
    this->S2.setTensileEnergy(this->S1.getTensileEnergy());
  }
};

#endif
