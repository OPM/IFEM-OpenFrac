// $Id$
//==============================================================================
//!
//! \file QuasiStaticSIM.h
//!
//! \date Sep 22 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Quasi-static solution driver for fracture elasticity simulators.
//!
//==============================================================================

#ifndef _QUASI_STATIC_SIM_H
#define _QUASI_STATIC_SIM_H

#include "NonLinSIM.h"


/*!
  \brief Nonlinear quasi-static solution driver for fracture simulators.
  \details This sub-class reimplements the lineSearch method which is
  taylored for the monolithic fracture elasticity problem.
*/

class QuasiStaticSIM : public NonLinSIM
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  QuasiStaticSIM(SIMbase& sim);
  //! \brief Empty destructor.
  virtual ~QuasiStaticSIM() {}

  //! \brief Prints out problem-specific data to the log stream.
  virtual void printProblem() const;

protected:
  //! \brief Performs line search to accelerate convergence.
  virtual bool lineSearch(TimeStep& param);
};

#endif
