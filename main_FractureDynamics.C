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
#include "SIMExplPhaseField.h"
#include "SIMFractureQstatic.h"
#ifdef IFEM_HAS_POROELASTIC
#include "SIMPoroElasticity.h"
#endif
#include "SIMCoupledSI.h"
#include "SIMSolverTS.h"
#include "FractureArgs.h"
#include "FractureElasticityVoigt.h"
#include "HHTSIM.h"
#include "GenAlphaSIM.h"
#include "NewmarkNLSIM.h"
#include "NonLinSIM.h"
#include "ASMstruct.h"
#include "ASMmxBase.h"


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

  //! \brief Overrides the stop time that was read from the input file.
  void setStopTime(double t) { Solver<T>::tp.stopTime = t; }

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
  \param[in] stopTime Stop time of the simulation (if non-negative)
  \param[in] context Input-file context for the time integrator
*/

template<class ElSolver, class PhaseSolver, class SIMFractureDynamics,
         template<class T1> class Solver=SIMSolver>
int runCombined (char* infile, double stopTime, const char* context)
{
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  ElSolver elastoSim;
  int rest = elastoSim.restartBasis(IFEM::getOptions().restartFile,
                                    IFEM::getOptions().restartStep);
  if (rest < 0) return -rest;

  ASMstruct::resetNumbering();
  if (!readModel(elastoSim,infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  PhaseSolver phaseSim(&elastoSim);
  if (!readModel(phaseSim,infile))
    return 1;

  phaseSim.opt.print(IFEM::cout) << std::endl;

  SIMFractureDynamics frac(elastoSim,phaseSim,infile);
  SIMDriver<SIMFractureDynamics,Solver> solver(frac,context);
  if (!readModel(solver,infile))
    return 1;

  if (stopTime >= 0.0)
    solver.setStopTime(stopTime);

  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
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
  \brief Determines whether the quasi-static semi-implicit driver is to be used.
*/

template<class ElSolver, class PhaseSolver,
         template<class T1, class T2> class Cpl,
         template<class T1> class Solver=SIMSolver>
int runSimulator6 (const FractureArgs& args, const char* context)
{
  if (args.integrator == 3 && args.coupling == 2)
  {
    typedef SIMFractureQstatic<ElSolver,PhaseSolver> Coupler;
    return runCombined<ElSolver,PhaseSolver,Coupler,Solver>(args.inpfile,
                                                            args.stopT,
                                                            context);
  }
  else if (args.integrator == 3 && args.coupling == 3)
  {
    typedef SIMFractureMiehe<ElSolver,PhaseSolver> Coupler;
    return runCombined<ElSolver,PhaseSolver,Coupler,Solver>(args.inpfile,
                                                            args.stopT,
                                                            context);
  }
  else
  {
    typedef SIMFracture<ElSolver,PhaseSolver,Cpl> Coupler;
    return runCombined<ElSolver,PhaseSolver,Coupler,Solver>(args.inpfile,
                                                            args.stopT,
                                                            context);
  }
}


/*!
  \brief Determines whether the explicit phase-field driver is to be used.
*/

template<class Dim, class ElSolver,
         template<class T1, class T2> class Cpl,
         template<class T1> class Solver=SIMSolver>
int runSimulator5 (const FractureArgs& args, const char* context)
{
  if (args.expPhase)
    return runSimulator6<ElSolver,SIMExplPhaseField,Cpl,Solver>(args,context);
  else
    return runSimulator6<ElSolver,SIMPhaseField<Dim>,Cpl,Solver>(args,context);
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
  typedef SIMDynElasticity<Dim,Integrator,ElSolver> SIMElastoDynamics;

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
         template<class T1, class T2> class Cpl>
int runSimulator3 (const FractureArgs& args)
{
  typedef SIMDynElasticity<Dim,Integrator,ElSolver> DynElSolver;

  const char* context = Integrator::inputContext;

  if (args.adap)
    return runSimulator5<Dim,DynElSolver,Cpl,SIMSolverTS>(args,context);

  return runSimulator5<Dim,DynElSolver,Cpl>(args,context);
}


/*!
  \brief Creates the combined fracture simulator and launches the simulation.
*/

template<class Dim, class Integrator,
         template<class T1, class T2> class Cpl>
int runSimulator2 (const FractureArgs& args)
{
  if (args.poroEl)
#ifdef IFEM_HAS_POROELASTIC
    return runSimulator3<Dim,Integrator,SIMPoroElasticity<Dim>,Cpl>(args);
#else
    return 99; // Built without the poroelastic coupling
#endif

  return runSimulator3<Dim,Integrator,SIMElasticityWrap<Dim>,Cpl>(args);
}


/*!
  \brief Selects the coupling driver to be used.
*/

template<class Dim, class Integrator>
int runSimulator1 (const FractureArgs& args)
{
  switch (args.coupling) {
  case 1:
  case 3:
    return runSimulator2<Dim,Integrator,SIMCoupled>(args);
  case 2:
    return runSimulator2<Dim,Integrator,SIMCoupledSI>(args);
  default: // No phase field coupling
    return runSimulator4<Dim,Integrator>(args,Integrator::inputContext);
  }
}


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

  FractureArgs args;
  Elastic::planeStrain = true;
  ASMmxBase::Type = ASMmxBase::NONE;
  IFEM::Init(argc,argv,"Fracture dynamics solver");

  for (int i = 1; i < argc; i++)
    if (argv[i] == args.inpfile || args.parseArg(argv[i]))
      ; // ignore the input file on the second pass
    else if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-mixed"))
      ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
    else if (!strcmp(argv[i],"-principal"))
      Elasticity::wantPrincipalStress = true;
    else if (!strncmp(argv[i],"-dbgEl",6) && i < argc-1)
      FractureElasticNorm::dbgElm = atoi(argv[++i]);
    else if (!strncmp(argv[i],"-stopT",6) && i < argc-1)
      args.stopT = atof(argv[++i]);
    else if (!args.inpfile)
      args.parseFile(argv[i],i);
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

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
