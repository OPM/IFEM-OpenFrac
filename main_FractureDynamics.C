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

#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMDynElasticity.h"
#include "SIMPhaseField.h"
#include "SIMFractureQstatic.h"
#include "SIMCoupledSI.h"
#include "SIMSolver.h"
#include "SIMSolverTS.h"
#include "HHTSIM.h"
#include "GenAlphaSIM.h"
#include "NewmarkNLSIM.h"
#include "NonLinSIM.h"
#include "ASMstruct.h"


//! \brief A struct collecting the command-line argument values.
struct FDargs
{
  //! Time integrator to use (0=linear quasi-static, no phase-field coupling,
  //! 1=linear Newmark, 2=linear generalized alpha, 3=nonlinear quasi-static,
  //! 4=nonlinear Hilber-Hughes-Taylor)
  char integrator;
  char coupling;  //!< Coupling flag (0: none, 1: staggered, 2: semi-implicit)
  bool adaptive;  //!< If \e true, use the time-slab adaptive solver
  char* inpfile;  //!< The input file to parse

  //! \brief Default constructor.
  FDargs() : inpfile(nullptr)
  { coupling = integrator = 1; adaptive = false; }
};


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
  //! \brief Empty destructor.
  virtual ~SIMDriver() {}

protected:
  using Solver<T>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),context))
    {
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        if (!strncasecmp(child->Value(),"stag",4))
          this->S1.parseStaggering(child);
        else
          this->SIMSolver<T>::parse(child);
    }
    else if (!strcasecmp(elem->Value(),"postprocessing"))
    {
      const TiXmlElement* child = elem->FirstChildElement("energyfile");
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
  \param[in] context Input-file context for the time integrator
*/

template<class ElSolver, class PhaseSolver, class SIMFractureDynamics,
         template<class T1> class Solver=SIMSolver>
int runCombined (char* infile, const char* context)
{
  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  ElSolver elastoSim;
  ASMstruct::resetNumbering();
  if (!elastoSim.read(infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  PhaseSolver phaseSim(&elastoSim);
  if (!phaseSim.read(infile))
    return 1;

  phaseSim.opt.print(IFEM::cout) << std::endl;

  SIMFractureDynamics frac(elastoSim,phaseSim,infile);
  SIMDriver<SIMFractureDynamics,Solver> solver(frac,context);
  if (!solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");
  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!frac.preprocess())
    return 2;

  // Initialize the linear solvers
  if (!elastoSim.initSystem(elastoSim.opt.solver,1,1,0,true) ||
      !phaseSim.initSystem(phaseSim.opt.solver))
    return 2;

  // Initialize the solution fields
  frac.init(TimeStep());

  if (elastoSim.opt.dumpHDF5(infile))
    solver.handleDataOutput(elastoSim.opt.hdf5,elastoSim.opt.saveInc);

  frac.setupDependencies();

  return solver.solveProblem(infile,"100. Starting the simulation",
                             phaseSim.getInitRefine() < 1);
}


/*!
  \brief Determines whether the quasi-static semi-implicit driver is to be used.
*/

template<class Dim, class Integrator,
         template<class T1, class T2> class Cpl,
         template<class T1> class Solver=SIMSolver>
int runSimulator3 (const FDargs& args)
{
  typedef SIMDynElasticity<Dim,Integrator> ElSolver;
  typedef SIMPhaseField<Dim>               PhaseSolver;

  const char* contx = Integrator::inputContext;
  if (args.integrator == 3 && args.coupling == 2)
  {
    typedef SIMFractureQstatic<ElSolver,PhaseSolver> Coupler;
    return runCombined<ElSolver,PhaseSolver,Coupler,Solver>(args.inpfile,contx);
  }
  else
  {
    typedef SIMFracture<ElSolver,PhaseSolver,Cpl> Coupler;
    return runCombined<ElSolver,PhaseSolver,Coupler,Solver>(args.inpfile,contx);
  }
}


/*!
  \brief Creates and launches a stand-alone elasticity simulator (no coupling).
  \param[in] infile The input file to parse
  \param[in] context Input-file context for the time integrator
*/

template<class Dim, class Integrator=NewmarkSIM>
int runStandAlone (char* infile, const char* context)
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

  // Initialize the solution fields
  elastoSim.init(TimeStep());

  if (elastoSim.opt.dumpHDF5(infile))
    solver.handleDataOutput(elastoSim.opt.hdf5,elastoSim.opt.saveInc);

  return solver.solveProblem(infile,"100. Starting the simulation");
}


/*!
  \brief Determines whether the adaptive simulation driver is to be used.
*/

template<class Dim, class Integrator, template<class T1, class T2> class Cpl>
int runSimulator2 (const FDargs& args)
{
  if (args.adaptive)
    return runSimulator3<Dim,Integrator,Cpl,SIMSolverTS>(args);

  return runSimulator3<Dim,Integrator,Cpl>(args);
}


/*!
  \brief Selects the coupling driver to be used.
*/

template<class Dim, class Integrator>
int runSimulator1 (const FDargs& args)
{
  switch (args.coupling) {
  case 1:
    return runSimulator2<Dim,Integrator,SIMCoupled>(args);
  case 2:
    return runSimulator2<Dim,Integrator,SIMCoupledSI>(args);
  default: // No phase field coupling
    return runStandAlone<Dim,Integrator>(args.inpfile,Integrator::inputContext);
  }
}


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
  \brief Selects the time integration driver to be used.
*/

template<class Dim>
int runSimulator (const FDargs& args)
{
  switch (args.integrator) {
  case 0:
    return runStandAlone<Dim,LinSIM>(args.inpfile,"staticsolver");
  case 1:
    return runSimulator1<Dim,NewmarkSIM>(args);
  case 2:
    return runSimulator1<Dim,GenAlphaSIM>(args);
  case 3:
    return runSimulator1<Dim,NonLinSIM>(args);
  case 4:
    return runSimulator1<Dim,HHTSIM>(args);
  case 5:
    return runSimulator1<Dim,NewmarkNLSIM>(args);
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

  FDargs args;
  bool twoD = false;

  IFEM::Init(argc,argv,"Fracture dynamics solver");

  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      twoD = SIMElasticity<SIM2D>::planeStrain = true;
    else if (!strcmp(argv[i],"-nocrack"))
      args.coupling = 0;
    else if (!strcmp(argv[i],"-semiimplicit"))
      args.coupling = 2;
    else if (!strcmp(argv[i],"-lstatic"))
      args.integrator = 0;
    else if (!strcmp(argv[i],"-GA"))
      args.integrator = 2;
    else if (!strcmp(argv[i],"-qstatic"))
      args.integrator = 3;
    else if (!strcmp(argv[i],"-HHT"))
      args.integrator = 4;
    else if (!strcmp(argv[i],"-oldHHT"))
      args.integrator = 5;
    else if (!strcmp(argv[i],"-principal"))
      Elasticity::wantPrincipalStress = true;
    else if (!strcmp(argv[i],"-dbgElm") && i < argc-1)
      FractureElasticNorm::dbgElm = atoi(argv[++i]);
    else if (!strncmp(argv[i],"-adap",5))
      args.adaptive = true;
    else if (!args.inpfile)
      args.inpfile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!args.inpfile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D] [-nGauss <n>]\n       "
              <<"[-nocrack|-semiimplicit] [-[l|q]static|-GA|-HHT] [-adaptive]\n"
              <<"       [-vtf <format> [-nviz <nviz>] [-nu <nu>] [-nv <nv]"
              <<" [-nw <nw>]] [-hdf5] [-principal]\n";
    return 0;
  }

  if (args.adaptive)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\nInput file: "<< args.inpfile;
  IFEM::getOptions().print(IFEM::cout) << std::endl;

  if (twoD)
    return runSimulator<SIM2D>(args);
  else
    return runSimulator<SIM3D>(args);
}
