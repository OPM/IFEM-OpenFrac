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
#include "SIMCoupledSI.h"
#include "SIMSolver.h"
#include "SIMSolverTS.h"
#include "HHTSIM.h"
#include "GenAlphaSIM.h"
#include "NewmarkNLSIM.h"
#include "NonLinSIM.h"
#include "ASMstruct.h"
#include "AppCommon.h"


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
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),context))
    {
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
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
*/

template<class Dim, class Integrator,
         template<class T1, class T2> class Cpl,
         template<class T1> class Solver=SIMSolver>
int runSimulator2 (char* infile)
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
  SIMDriver<SIMFractureDynamics,Solver> solver(frac,Integrator::inputContext);
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
int runSimulator3 (char* infile, const char* context)
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


/*!
  \brief Creates the combined fracture simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] timeslabs Use time-slab adaptive solver
*/

template<class Dim, class Integrator, template<class T1, class T2> class Cpl>
int runSolver (char* infile, bool timeslabs)
{
  if (timeslabs)
    return runSimulator2<Dim,Integrator,Cpl,SIMSolverTS>(infile);

  return runSimulator2<Dim,Integrator,Cpl>(infile);
}


/*!
  \brief Creates the combined fracture simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] coupling Coupling flag (0: none, 1: staggered, 2: semi-implicit)
  \param[in] timeslabs Use time-slab adaptive solver
*/

template<class Dim, class Integrator=NewmarkSIM>
int runSimulator1 (char* infile, char coupling, bool timeslabs)
{
  if (coupling == 1)
    return runSolver<Dim,Integrator,SIMCoupled>(infile,timeslabs);
  else if (coupling == 2)
    return runSolver<Dim,Integrator,SIMCoupledSI>(infile,timeslabs);
  else // No phase field coupling
    return runSimulator3<Dim,Integrator>(infile,Integrator::inputContext);
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
  \brief Creates the combined fracture simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] integrator The time integrator to use (0=linear quasi-static,
             no phase-field coupling, 1=linear Newmark, 2=Generalized alpha,
             3=nonlinear quasi-static, 4=nonlinear Hilber-Hughes-Taylor)
  \param[in] coupling Coupling flag (0: none, 1: staggered, 2: semi-implicit)
  \param[in] timeslabs Use time-slab adaptive solver
*/

template<class Dim>
int runSimulator (char* infile, char integrator, char coupling, bool timeslabs)
{
  switch (integrator) {
  case 0:
    return runSimulator3<Dim,LinSIM>(infile,"staticsolver");
  case 1:
    return runSimulator1<Dim>(infile,coupling,timeslabs);
  case 2:
    return runSimulator1<Dim,GenAlphaSIM>(infile,coupling,timeslabs);
  case 3:
    return runSimulator1<Dim,NonLinSIM>(infile,coupling,timeslabs);
  case 4:
    return runSimulator1<Dim,HHTSIM>(infile,coupling,timeslabs);
  case 5:
    return runSimulator1<Dim,NewmarkNLSIM>(infile,coupling,timeslabs);
  default:
    std::cerr <<" *** Invalid time integrator "<< integrator << std::endl;
    return 99;
  }
}


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
    else if (!strcmp(argv[i],"-lstatic"))
      integrator = 0;
    else if (!strcmp(argv[i],"-GA"))
      integrator = 2;
    else if (!strcmp(argv[i],"-qstatic"))
      integrator = 3;
    else if (!strcmp(argv[i],"-HHT"))
      integrator = 4;
    else if (!strcmp(argv[i],"-oldHHT"))
      integrator = 5;
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
              <<"       [-lag|-spec|-LR] [-2D] [-nGauss <n>]\n       "
              <<"[-nocrack|-semiimplicit] [-[l|q]static|-GA|-HHT] [-adaptive]\n"
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
    return runSimulator<SIM2D>(infile,integrator,coupling,adaptive);
  else
    return runSimulator<SIM3D>(infile,integrator,coupling,adaptive);
}
