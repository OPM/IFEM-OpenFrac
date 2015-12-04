// $Id$
//==============================================================================
//!
//! \file main_FractureDynamics.C
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving
//!
//! \brief Main program for an isogeometric fracture-dynamics solver.
//!
//==============================================================================

#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMPhaseField.h"
#include "SIMElasticityWrap.h"
#include "SIMCoupled.h"
#include "SIMSolver.h"
#include "ASMstruct.h"
#include "AppCommon.h"


template<class Dim, class Integrand> int runSimulator2 (char* infile)
{
  typedef SIMPhaseField<Dim,Integrand>                SIMCrackField;
  typedef SIMElasticityWrap<Dim>                      SIMElastoDynamics;
  typedef SIMCoupled<SIMElastoDynamics,SIMCrackField> SIMFractureDynamics;

  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  SIMCrackField phaseSim;
  if (!phaseSim.read(infile))
    return 1;

  phaseSim.opt.print(IFEM::cout) << std::endl;

  SIMElastoDynamics elastoSim;
  ASMstruct::resetNumbering();
  if (!elastoSim.read(infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  SIMFractureDynamics            frac(elastoSim,phaseSim);
  SIMSolver<SIMFractureDynamics> solver(frac);
  if (!solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");
  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!phaseSim.preprocess() || !elastoSim.preprocess())
    return 2;

  // Initialize the linear solvers
  phaseSim.initSystem(phaseSim.opt.solver);
  elastoSim.initSystem(elastoSim.opt.solver);
  elastoSim.initSol();

  // Time-step loop
  phaseSim.init(TimeStep());
  phaseSim.setInitialConditions();

  DataExporter* exporter = NULL;
  if (phaseSim.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(frac,solver,phaseSim.opt.hdf5,false,1,1);

  phaseSim.registerDependency(&elastoSim,"tensile",1);
  // This is defined on integration point and not on control points.
  // It is a global vector across all patches on the process.
  // Use an explicit call instead of normal couplings for this.
  phaseSim.setTensileEnergy(elastoSim.getTensileEnergy()->data());

  int res = solver.solveProblem(infile,exporter,"100. Starting the simulation");

  delete exporter;
  return res;
}


template<class Dim> int runSimulator1 (char* infile, bool fourth)
{
  if (fourth)
    return runSimulator2<Dim,CahnHilliard4>(infile);
  else
    return runSimulator2<Dim,CahnHilliard>(infile);
}


int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  int  i;
  char* infile = 0;
  bool twoD = false;
  bool fourth = false;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      twoD = true;
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
              <<" [-lag|-spec|-LR] [-2D] [-nGauss <n>]"
              <<"\n       [-vtf <format> [-nviz <nviz>]"
              <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n"<< std::endl;
    return 0;
  }

  IFEM::cout <<"\n >>> IFEM Fracture dynamics solver <<<"
             <<"\n =====================================\n"
             <<"\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
  IFEM::cout <<"\n\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (twoD)
    return runSimulator1<SIM2D>(infile,fourth);
  else
    return runSimulator1<SIM3D>(infile,fourth);
}
