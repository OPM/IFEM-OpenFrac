// $Id$
//==============================================================================
//!
//! \file FractureLaunch.h
//!
//! \date Nov 27 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Launching the simulators of the fracture-dynamics solver.
//! \details The launch of each time integrator is instantiated for each
//! dimension in a translation unit of its own,
//! FractureLaunch\<Integrator\>\<Dim\>D.C, as each one is a product of the
//! couplings, elasticity models, phase fields and adaptivity.
//!
//==============================================================================

#ifndef _FRACTURE_LAUNCH_H
#define _FRACTURE_LAUNCH_H

#include "IFEM.h"
#include "Profiler.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMDynElasticity.h"
#include "SIMPhaseField.h"
#include "SIMExplPhaseField.h"
#include "SIMFractureQstatic.h"
#ifdef IFEM_HAS_POROELASTIC
#include "SIMPoroElasticity.h"
#else
#include "SIMElasticityWrap.h"
#endif
#include "SIMCoupledSI.h"
#include "SIMSolverTS.h"
#include "ElasticityUtils.h"
#include "FractureArgs.h"
#include "FractureElasticityVoigt.h"
#include "HHTSIM.h"
#include "GenAlphaSIM.h"
#include "NewmarkNLSIM.h"
#include "NonLinSIM.h"
#include "ASMstruct.h"
#include "ASMmxBase.h"

#include <type_traits>


/*!
  \brief Dynamic simulation driver.

  \details Only the parse method is reimplemented here to handle that the
  time stepping parameters may be located within the specified context.
*/

template<class T, template<class S1> class Solver>
class SIMDriver : public Solver<T>
{
public:
  //! \brief The constructor initializes the reference to the actual solver.
  SIMDriver(T& s, const char* c = nullptr) : Solver<T>(s), context(c) {}

  //! \brief Overrides the stop time that was read from the input file.
  void setStopTime(double t) { Solver<T>::tp.stopTime = t; }

protected:
  using Solver<T>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem)
  {
    if (!strcasecmp(elem->Value(),context))
    {
      const tinyxml2::XMLElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        if (!strncasecmp(child->Value(),"stag",4))
          this->S1.parseStaggering(child);
        else
          this->SIMSolver<T>::parse(child);
    }
    else if (!strcasecmp(elem->Value(),"postprocessing"))
    {
      const tinyxml2::XMLElement* child = elem->FirstChildElement("energyfile");
      if (child && child->FirstChild())
        this->S1.setEnergyFile(child->FirstChild()->Value());
    }

    return this->Solver<T>::parse(elem);
  }

private:
  const char* context; //!< XML-tag to search for time-stepping input within
};


/*!
  \brief Creates the combined fracture simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] stopTime Stop time of the simulation (if non-negative)
  \param[in] context Input-file context for the time integrator
*/

template<class ElSolver, class PhaseSolver, class FractureSim,
         template<class T1> class Solver=SIMSolver>
int runCombined (char* infile, double stopTime, const char* context)
{
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  ElSolver elastoSim;
  if (int rest = elastoSim.restartBasis(IFEM::getOptions().restartFile,
                                        IFEM::getOptions().restartStep);
      rest < 0) return -rest;

  ASMstruct::resetNumbering();
  if (!readModel(elastoSim,infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  PhaseSolver phaseSim(&elastoSim);
  if (!readModel(phaseSim,infile))
    return 1;

  phaseSim.opt.print(IFEM::cout) << std::endl;

  FractureSim frac(elastoSim,phaseSim,infile);
  SIMDriver<FractureSim,Solver> solver(frac,context);
  if (!readModel(solver,infile))
    return 1;

  if (stopTime >= 0.0)
    solver.setStopTime(stopTime);

  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  elastoSim.initPrm();
  if (!frac.preprocess())
    return 2;

  // Initialize the linear solvers
  if (!elastoSim.initSystem(elastoSim.opt.solver,1,1,1,true) ||
      !phaseSim.initSystem(phaseSim.opt.solver,1,1,1))
    return 2;

  // Initialize the solution fields
  if (!frac.init(solver.getTimePrm()))
    return 2;

  if (solver.restart(elastoSim.opt.restartFile,elastoSim.opt.restartStep) < 0)
    return 2;

  if (elastoSim.opt.dumpHDF5(infile))
    solver.handleDataOutput(elastoSim.opt.hdf5,elastoSim.getProcessAdm(),
                            elastoSim.opt.saveInc,elastoSim.opt.restartInc);

  frac.setupDependencies();

  return solver.solveProblem(infile,"100. Starting the simulation");
}


/*!
  \brief Creates and launches a stand-alone elasticity simulator (no coupling).
  \param[in] infile The input file to parse
  \param[in] stopTime Stop time of the simulation (if non-negative)
  \param[in] context Input-file context for the time integrator
*/

template<class Dim, class Integrator, class ElSolver>
int runStandAlone (char* infile, double stopTime, const char* context)
{
  using SIMElastoDynamics = SIMDynElasticity<Dim,Integrator,ElSolver>;

  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  SIMElastoDynamics elastoSim;
  if (!readModel(elastoSim,infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  SIMDriver<SIMElastoDynamics,SIMSolver> solver(elastoSim,context);
  if (!readModel(solver,infile))
    return 1;

  if (stopTime >= 0.0)
    solver.setStopTime(stopTime);

  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!elastoSim.preprocess())
    return 2;

  // Initialize the linear solvers
  if (!elastoSim.initSystem(elastoSim.opt.solver,1,1,0,true))
    return 2;

  // Initialize the solution fields
  if (!elastoSim.init(TimeStep()))
    return 2;

  if (solver.restart(elastoSim.opt.restartFile,elastoSim.opt.restartStep) < 0)
    return 2;

  if (elastoSim.opt.dumpHDF5(infile))
    solver.handleDataOutput(elastoSim.opt.hdf5,elastoSim.getProcessAdm(),
                            elastoSim.opt.saveInc,elastoSim.opt.restartInc);

  return solver.solveProblem(infile,"100. Starting the simulation");
}


/*!
  \brief Determines whether the explicit phase-field driver is to be used.
*/

template<class Dim, class ElSolver,
         template<class T1, class T2> class Cpl,
         template<class T1, class T2,
                  template<class T3, class T4> class Coupling> class FractureSim,
         template<class T1> class Solver=SIMSolver>
int runSimulator5 (const FractureArgs& args, const char* context)
{
  if (args.expPhase) {
    using FracSim = FractureSim<ElSolver,SIMExplPhaseField,Cpl>;
    return runCombined<ElSolver,SIMExplPhaseField,FracSim,Solver>(args.inpfile,
                                                                  args.stopT,
                                                                  context);

  } else {
    using FracSim = FractureSim<ElSolver,SIMPhaseField<Dim>,Cpl>;
    return runCombined<ElSolver,SIMPhaseField<Dim>,FracSim,Solver>(args.inpfile,
                                                                   args.stopT,
                                                                   context);
  }
}


/*!
  \brief Determines whether the poroelastic simulation driver is to be used.
*/

template<class Dim, class Integrator>
int runSimulator4 (const FractureArgs& args,
                   const char* context = "newmarksolver")
{
  if (args.poroEl)
#ifdef IFEM_HAS_POROELASTIC
    return runStandAlone<Dim,Integrator,SIMPoroElasticity<Dim>>(args.inpfile,
                                                                args.stopT,
                                                                context);
#else
    return 99; // Built without the poroelastic coupling
#endif

  return runStandAlone<Dim,Integrator,SIMElasticityWrap<Dim>>(args.inpfile,
                                                              args.stopT,
                                                              context);
}


/*!
  \brief Determines whether the adaptive simulation driver is to be used.
*/

template<class Dim, class Integrator, class ElSolver,
         template<class T1, class T2> class Cpl,
         template<class T1, class T2,
                  template<class T3, class T4> class Coupling> class FractureSim>
int runSimulator3 (const FractureArgs& args)
{
  using DynElSolver = SIMDynElasticity<Dim,Integrator,ElSolver>;

  const char* context = Integrator::inputContext;

  if (args.adap)
    return runSimulator5<Dim,DynElSolver,Cpl,FractureSim,SIMSolverTS>(args,context);

  return runSimulator5<Dim,DynElSolver,Cpl,FractureSim>(args,context);
}


/*!
  \brief Creates the combined fracture simulator and launches the simulation.
*/

template<class Dim, class Integrator,
         template<class T1, class T2> class Cpl,
         template<class T1, class T2,
                  template<class T3, class T4> class Coupling> class FractureSim>
int runSimulator2 (const FractureArgs& args)
{
  if (args.poroEl)
#ifdef IFEM_HAS_POROELASTIC
    return runSimulator3<Dim,Integrator,SIMPoroElasticity<Dim>,Cpl,FractureSim>(args);
#else
    return 99; // Built without the poroelastic coupling
#endif

  return runSimulator3<Dim,Integrator,SIMElasticityWrap<Dim>,Cpl,FractureSim>(args);
}


//! \brief Tells whether a fracture simulator is the generic one.
template<template<class T1, class T2,
                  template<class T3, class T4> class Coupling> class FractureSim>
struct IsSIMFracture : std::false_type {};

//! \brief Specialization for the generic fracture simulator.
template<>
struct IsSIMFracture<SIMFracture> : std::true_type {};


/*!
  \brief Selects the coupling driver to be used.
  \details The nonlinear quasi-static simulation has fracture simulators of
  its own for the couplings 2 and 3, see runSimulator, so the generic
  simulator is not instantiated for the semi-implicit coupling there.
*/

template<class Dim, class Integrator,
         template<class T1, class T2,
                  template<class T3, class T4> class Coupling> class FractureSim=SIMFracture>
int runSimulator1 (const FractureArgs& args)
{
  switch (args.coupling) {
  case 1:
  case 3:
    return runSimulator2<Dim,Integrator,SIMCoupled,FractureSim>(args);
  case 2:
    if constexpr (std::is_same_v<Integrator,NonLinSIM> &&
                  IsSIMFracture<FractureSim>::value)
      return 99; // Never reached, coupling 2 is run by SIMFractureQstatic
    else
      return runSimulator2<Dim,Integrator,SIMCoupledSI,FractureSim>(args);
  default: // No phase field coupling
    return runSimulator4<Dim,Integrator>(args,Integrator::inputContext);
  }
}


//! \brief Helper 3-param Qstatic template.
template<class T1, class T2, template<class T3, class T4> class Dummy>
using QstaticFracture = SIMFractureQstatic<T1,T2>;

//! \brief Helper 3-param Miehe template.
template<class T1, class T2, template<class T3, class T4> class Dummy>
using MieheFracture = SIMFractureMiehe<T1,T2>;


/*!
  \brief Launches the simulation with a time integration driver.
  \details The nonlinear quasi-static simulation has fracture simulators of
  its own for the couplings 2 and 3.
*/

template<class Dim, class Integrator>
int runIntegrator (const FractureArgs& args);

#endif
