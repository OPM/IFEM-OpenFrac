// $Id$
//==============================================================================
//!
//! \file main_FractureDynamics.C
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving
//!
//! \brief Main program for the isogeometric fracture-dynamics solver.
//!
//==============================================================================

#include "FractureLaunch.h"


/*!
  \brief Selects the time integration driver to be used.
*/

template<class Dim>
int runSimulator (const FractureArgs& args)
{
  switch (args.integrator) {
  case 0:
    if (args.coupling > 0)
    {
      std::cerr <<" *** The linear static option is for non-cracking material"
                <<" only."<< std::endl;;
      return 98;
    }
    return runSimulator4<Dim,LinSIM>(args,"staticsolver");
  case 1:
    return runIntegrator<Dim,NewmarkSIM>(args);
  case 2:
    return runIntegrator<Dim,GenAlphaSIM>(args);
  case 3:
    return runIntegrator<Dim,NonLinSIM>(args);
  case 4:
    return runIntegrator<Dim,HHTSIM>(args);
  case 5:
    return runIntegrator<Dim,NewmarkNLSIM>(args);
  default:
    std::cerr <<" *** Invalid time integrator "<< args.integrator << std::endl;
    return 99;
  }
}


/*!
  \brief Main program for the isogeometric fracture elasticity solver.
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);

  IFEM::Init(argc,argv,"Fracture dynamics solver");

  ASMmxBase::Type = ASMmxBase::NONE;
  Elastic::planeStrain = true;
  FractureArgs args;

  // Lambda function for parsing the options not stored in FractureArgs.
  auto external_args = [argc,argv,&args](int& i) -> bool
  {
    if (!strcmp(argv[i],"-mixed"))
      ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
    else if (!strcmp(argv[i],"-principal"))
      Elasticity::wantPrincipalStress = true;
    else if (!strncmp(argv[i],"-dbgEl",6) && i < argc-1)
      FractureElasticNorm::dbgElm = atoi(argv[++i]);
    else if (!strncmp(argv[i],"-stopT",6) && i < argc-1)
      args.stopT = atof(argv[++i]);
    else
      return false;
    return true;
  };

  args.inpfile = args.preParse(argc,argv,external_args);
  if (!args.inpfile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D] [-mixed] [-nGauss <n>]\n"
              <<"       [-nocrack|-explcrack|-semiimplicit]"
              <<" [-[l|q]static|-GA|-HHT] [-poro] [-adaptive]\n"
              <<"       [-vtf <format> [-nviz <nviz>] [-nu <nu>] [-nv <nv]"
              <<" [-nw <nw>]]\n       [-hdf5] [-principal] [-stopTime <t>]\n";
    return 0;
  }

  if (args.adap)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\nInput file: "<< args.inpfile;
  IFEM::getOptions().print(IFEM::cout);
  if (args.stopT >= 0.0)
    IFEM::cout <<"\nSimulation stop time: "<< args.stopT;
  IFEM::cout << std::endl;

  if (args.dim == 2)
    return runSimulator<SIM2D>(args);
  else if (args.dim == 3)
    return runSimulator<SIM3D>(args);

  return 1; // No 1D solution
}
