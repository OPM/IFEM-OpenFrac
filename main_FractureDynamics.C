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
#include "SIMDynElasticity.h"
#include "SIMPhaseField.h"
#include "SIMFractureDynamics.h"
#include "SIMSolver.h"
#include "SIMDriver.h"
#include "ASMstruct.h"
#include "AppCommon.h"


/*!
  \brief Creates the combined fracture simulator and launches the simulation.
  \param[in] infile The input file to parse
*/

template<class Dim, class Integrator,
         template<class T1, class T2> class Cpl,
         template<class T1> class Solver=SIMSolver>
int runCplSimulator (char* infile)
{
  typedef SIMDynElasticity<Dim,Integrator> SIMElastoDynamics;
  typedef SIMPhaseField<Dim>               SIMCrackField;

  typedef SIMFracture<SIMElastoDynamics,SIMCrackField,Cpl> SIMFractureDynamics;

  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  SIMElastoDynamics elastoSim;
  ASMstruct::resetNumbering();
  if (!elastoSim.read(infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  SIMCrackField phaseSim(&elastoSim);
  if (!phaseSim.read(infile))
    return 1;

  phaseSim.opt.print(IFEM::cout) << std::endl;

  SIMFractureDynamics frac(elastoSim,phaseSim,infile);
  SIMDriver<SIMFractureDynamics,Solver> solver(frac,"newmarksolver");
  if (!solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");
  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!frac.preprocess())
    return 2;

  // Initialize the linear solvers
  if (!elastoSim.initSystem(elastoSim.opt.solver) ||
      !phaseSim.initSystem(phaseSim.opt.solver,1,1,false))
    return 2;

  // Time-step loop
  frac.init(TimeStep());

  DataExporter* exporter = nullptr;
  if (elastoSim.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(frac,solver,elastoSim.opt.hdf5,false,1,1);

  frac.setupDependencies();

  int res = solver.solveProblem(infile,exporter,"100. Starting the simulation",
                                phaseSim.getInitRefine() < 1);

  delete exporter;
  return res;
}


/*!
  \brief Creates and launches a stand-alone elasticity simulator (no coupling).
  \param[in] infile The input file to parse
  \param[in] context Input-file context for the time integrator
*/

template<class Dim, class Integrator=NewmarkSIM>
int runSimulator (char* infile, const char* context = "newmarksolver")
{
  typedef SIMDynElasticity<Dim,Integrator> SIMElastoDynamics;

  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  SIMElastoDynamics elastoSim;
  if (!elastoSim.read(infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  SIMDriver<SIMElastoDynamics,SIMSolver> solver(elastoSim,context);
  if (!solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");
  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!elastoSim.preprocess())
    return 2;

  // Initialize the linear solvers
  if (!elastoSim.initSystem(elastoSim.opt.solver))
    return 2;

  // Time-step loop
  elastoSim.init(TimeStep());

  DataExporter* exporter = nullptr;
  if (elastoSim.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(elastoSim,solver,elastoSim.opt.hdf5,
                                     false,1,1);

  int res = solver.solveProblem(infile,exporter,"100. Starting the simulation");

  delete exporter;
  return res;
}


#include "runTimeInt.h"


/*!
  \brief Main program for NURBS-based fracture elasticity solver.
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);

  int  i;
  char* infile = 0;
  char coupling = 1;
  char integrator = 1;
  bool twoD = false;
  bool adaptive = false;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      twoD = SIMElasticity<SIM2D>::planeStrain = true;
    else if (!strcmp(argv[i],"-nocrack"))
      coupling = 0;
    else if (!strcmp(argv[i],"-semiimplicit"))
      coupling = 2;
    else if (!strcmp(argv[i],"-static"))
      integrator = 0;
    else if (!strcmp(argv[i],"-GA"))
      integrator = 2;
    else if (!strcmp(argv[i],"-principal"))
      Elasticity::wantPrincipalStress = true;
    else if (!strcmp(argv[i],"-dbgElm") && i < argc-1)
      FractureElasticNorm::dbgElm = atoi(argv[++i]);
    else if (!strncmp(argv[i],"-adap",5))
      adaptive = true;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D] [-nGauss <n>]\n"
              <<"       [-nocrack|-semiimplicit] [-static|-GA] [-adaptive]\n"
              <<"       [-vtf <format> [-nviz <nviz>] [-nu <nu>] [-nv <nv]"
              <<" [-nw <nw>]] [-hdf5] [-principal]\n"<< std::endl;
    return 0;
  }

  if (adaptive)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\n >>> IFEM Fracture dynamics solver <<<"
             <<"\n =====================================\n"
             <<"\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
  IFEM::cout <<"\n\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (twoD)
    return runTimeInt<SIM2D>(infile,integrator,coupling,adaptive);
  else
    return runTimeInt<SIM3D>(infile,integrator,coupling,adaptive);
}
