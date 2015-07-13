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
#include "SIMCH.h"
#include "Utilities.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "AppCommon.h"
#include "TimeStep.h"


template<class Dim>
int runSimulator(char* infile)
{
  SIMCH<Dim> ch;

  int res = ConfigureSIM(ch, infile, false);

  if (res)
    return res;

  // HDF5 output
  DataExporter* exporter=NULL;

  SIMSolver<SIMCH<Dim>> solver(ch);
  if (ch.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(ch, solver, ch.opt.hdf5,
                                     false, 1, 1);
  res = solver.solveProblem(infile, exporter);

  delete exporter;
  return res;
}


int main(int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  int  i;
  char ndim = 3;
  char* infile = 0;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
    else if (!strcmp(argv[i],"-1D"))
      ndim = 1;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
              <<" [-free] [-lag|-spec|-LR] [-1D|-2D] [-nGauss <n>]"
              <<"\n       [-vtf <format> [-nviz <nviz>]"
              <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n"
              <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
              <<"       [-ignore <p1> <p2> ...] [-fixDup]" << std::endl;
    return 0;
  }

  IFEM::cout <<"\n >>> IFEM Cahn-Hilliard equation solver <<<"
             <<"\n ====================================\n"
             <<"\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
  IFEM::cout <<"\n\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (ndim == 3)
    return runSimulator<SIM3D>(infile);
  else if (ndim == 2)
    return runSimulator<SIM2D>(infile);
  else
    return runSimulator<SIM1D>(infile);

  return 1;
}
