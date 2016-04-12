// $Id$
//==============================================================================
//!
//! \file runTimeInt.h
//!
//! \date Dec 1 2015
//!
//! \author Knut Morten Okstad
//!
//! \brief Coupled simulation drivers with time evolution.
//!
//==============================================================================

#ifndef _RUN_TIME_INT_H
#define _RUN_TIME_INT_H

#include "SIMCoupledSI.h"
#include "SIMSolverTS.h"
#include "GenAlphaSIM.h"
#include "NonLinSIM.h"


/*!
  \brief Linear quasi-static solution driver.
*/

class LinSIM : public NonLinSIM
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  LinSIM(SIMbase& sim) : NonLinSIM(sim,NonLinSIM::NONE) { fromIni = false; }
  //! \brief Empty destructor.
  virtual ~LinSIM() {}
};


/*!
  \brief Creates an adaptive simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] adaptive If \e true, use time-slab adaptive solver
*/

template<class Dim, class Integrator, template<class T1, class T2> class Cpl>
int runAdapSolver (char* infile, bool adaptive)
{
  if (adaptive)
    return runCplSimulator<Dim,Integrator,Cpl,SIMSolverTS>(infile);

  return runCplSimulator<Dim,Integrator,Cpl>(infile);
}


/*!
  \brief Creates a coupled fracture simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] coupling Coupling flag (0: none, 1: staggered, 2: semi-implicit)
  \param[in] adaptive If \e true, use time-slab adaptive solver
*/

template<class Dim, class Integrator=NewmarkSIM>
int runCplSolver (char* infile, char coupling, bool adaptive)
{
  if (coupling == 1)
    return runAdapSolver<Dim,Integrator,SIMCoupled>(infile,adaptive);
  else if (coupling == 2)
    return runAdapSolver<Dim,Integrator,SIMCoupledSI>(infile,adaptive);
  else
    return runSimulator<Dim,Integrator>(infile); // Single field, no coupling
}


/*!
  \brief Creates a simulator with time evolution and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] integrator The time integrator to use (0=linear quasi-static,
             no phase-field coupling, 1=linear Newmark, 2=Generalized alpha)
  \param[in] coupling Coupling flag (0: none, 1: staggered, 2: semi-implicit)
  \param[in] adaptive If \e true, use time-slab adaptive solver
*/

template<class Dim>
int runTimeInt (char* infile, char integrator, char coupling, bool adaptive)
{
  if (integrator == 2)
    return runCplSolver<Dim,GenAlphaSIM>(infile,coupling,adaptive);
  else if (integrator > 0)
    return runCplSolver<Dim>(infile,coupling,adaptive);
  else
    return runSimulator<Dim,LinSIM>(infile,"staticsolver"); // Single field
}

#endif
