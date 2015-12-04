// $Id$
//==============================================================================
//!
//! \file main_CahnHilliard.C
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving
//!
//! \brief Main program for an isogeometric solver for Cahn-Hilliard.
//!
//==============================================================================

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMPhaseField.h"
#include "SIMSolver.h"
#include "AppCommon.h"


template<class Dim, class Integrand> int runSimulator2 (char* infile)
{
  typedef SIMPhaseField<Dim,Integrand> PhaseFieldDriver;

  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  PhaseFieldDriver phaseSim;
  if (!phaseSim.read(infile))
    return 1;

  phaseSim.opt.print(IFEM::cout) << std::endl;

  SIMSolver<PhaseFieldDriver> solver(phaseSim);
  if (!solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");
  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!phaseSim.preprocess())
    return 2;

  // Initialize the linear solvers
  if (!phaseSim.initSystem(phaseSim.opt.solver,1,1,false))
    return 2;

  // Time-step loop
  phaseSim.init(TimeStep());
  phaseSim.setInitialConditions();

  DataExporter* exporter = NULL;
  if (phaseSim.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(phaseSim,solver,phaseSim.opt.hdf5,
                                     false,1,1);

  int res = solver.solveProblem(infile,exporter,"100. Starting the simulation");

  delete exporter;
  return res;
}


template<class Dim> int runSimulator1 (char* infile, bool fourth)
{
  if (fourth)
    return runSimulator2<Dim,CahnHilliard4>(infile);

  return runSimulator2<Dim,CahnHilliard>(infile);
}


int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  int i, ndim = 3;
  char* infile = 0;
  bool fourth = false;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
    else if (!strcmp(argv[i],"-1D"))
      ndim = 1;
    else if (!strcmp(argv[i],"-fourth"))
      fourth = true;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
              <<" [-lag|-spec|-LR] [-1D|-2D] [-nGauss <n>]"
              <<"\n       [-vtf <format> [-nviz <nviz>]"
              <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]"<< std::endl;
    return 0;
  }

  IFEM::cout <<"\n >>> IFEM Cahn-Hilliard equation solver <<<"
             <<"\n ==========================================\n"
             <<"\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
  IFEM::cout <<"\n\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (ndim == 3)
    return runSimulator1<SIM3D>(infile,fourth);
  else if (ndim == 2)
    return runSimulator1<SIM2D>(infile,fourth);
  else
    return runSimulator1<SIM1D>(infile,fourth);
}
