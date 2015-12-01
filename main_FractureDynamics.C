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
#include "SIMCoupled.h"
#include "SIMSolver.h"
#include "GenAlphaSIM.h"
#include "NonLinSIM.h"
#include "ASMstruct.h"
#include "AppCommon.h"


/*!
  \brief Dynamic simulation driver.

  \details Only the parse method is reimplemented here to handle that the
  time stepping parameters may be located within the specified context.
*/

template<class T> class SIMDriver : public SIMSolver<T>
{
public:
  //! \brief The constructor initializes the reference to the actual solver.
  SIMDriver(T& s, const char* c = nullptr) : SIMSolver<T>(s), context(c) {}
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

    return this->SIMSolver<T>::parse(elem);
  }

private:
  const char* context; //!< XML-tag to search for time-stepping input within
};


/*!
  \brief Creates the combined fracture simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] pfOrder Order of the crack phase field (2 or 4)
*/

template<class Dim, class Integrator=NewmarkSIM>
int runSimulator2 (char* infile, char pfOrder)
{
  typedef SIMDynElasticity<Dim,Integrator>            SIMElastoDynamics;
  typedef SIMPhaseField<Dim>                          SIMCrackField;
  typedef SIMCoupled<SIMElastoDynamics,SIMCrackField> SIMFractureDynamics;

  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  SIMElastoDynamics elastoSim;
  ASMstruct::resetNumbering();
  if (!elastoSim.read(infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  SIMCrackField phaseSim(pfOrder);
  if (!phaseSim.read(infile))
    return 1;

  phaseSim.opt.print(IFEM::cout) << std::endl;

  SIMFractureDynamics            frac(elastoSim,phaseSim);
  SIMDriver<SIMFractureDynamics> solver(frac,"newmarksolver");
  if (!solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");
  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!elastoSim.preprocess() || !phaseSim.preprocess())
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

  elastoSim.registerDependency(&phaseSim,"phasefield",1);
  // The tensile energy is defined on integration points and not control points.
  // It is a global buffer array across all patches in the model.
  // Use an explicit call instead of normal couplings for this.
  phaseSim.setTensileEnergy(elastoSim.getTensileEnergy());

  int res = solver.solveProblem(infile,exporter,
                                "100. Starting the simulation",false);

  delete exporter;
  return res;
}


/*!
  \brief Creates and launches a stand-alone elasticity simulator (no coupling).
  \param[in] infile The input file to parse
  \param[in] context Input-file context for the time integrator
*/

template<class Dim, class Integrator=NewmarkSIM>
int runSimulator3 (char* infile, const char* context = "newmarksolver")
{
  typedef SIMDynElasticity<Dim,Integrator> SIMElastoDynamics;

  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  SIMElastoDynamics elastoSim;
  if (!elastoSim.read(infile))
    return 1;

  elastoSim.opt.print(IFEM::cout) << std::endl;

  SIMDriver<SIMElastoDynamics> solver(elastoSim,context);
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
  \param[in] pfOrder Order of the crack phase field, 0=(no crack), 2 or 4
*/

template<class Dim, class Integrator=NewmarkSIM>
int runSimulator1 (char* infile, char phaseFieldOrder)
{
  if (phaseFieldOrder > 0)
    return runSimulator2<Dim,Integrator>(infile,phaseFieldOrder);
  else
    return runSimulator3<Dim,Integrator>(infile); // No phase field coupling
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
             no phase-field coupling, 1=linear Newmark, 2=Generalizes alpha)
  \param[in] pfOrder Order of the crack phase field, 0=(no crack), 2 or 4
*/

template<class Dim>
int runSimulator (char* infile, char integrator, char phaseFieldOrder)
{
  if (integrator == 2)
    return runSimulator1<Dim,GenAlphaSIM>(infile,phaseFieldOrder);
  else if (integrator > 0)
    return runSimulator1<Dim>(infile,phaseFieldOrder);
  else
    return runSimulator3<Dim,LinSIM>(infile,"staticsolver");
}


/*!
  \brief Main program for NURBS-based fracture elasticity solver.
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);

  int  i;
  char* infile = 0;
  char pfOrder = 2;
  char integrator = 1;
  bool twoD = false;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      twoD = SIMElasticity<SIM2D>::planeStrain = true;
    else if (!strcmp(argv[i],"-fourth"))
      pfOrder = 4;
    else if (!strcmp(argv[i],"-nocrack"))
      pfOrder = 0;
    else if (!strcmp(argv[i],"-static"))
      integrator = 0;
    else if (!strcmp(argv[i],"-GA"))
      integrator = 2;
    else if (!strcmp(argv[i],"-principal"))
      Elasticity::wantPrincipalStress = true;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D] [-nGauss <n>]\n"
              <<"       [-fourth|-nocrack] [-GA|-static]\n"
              <<"       [-vtf <format> [-nviz <nviz>]"
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
    return runSimulator<SIM2D>(infile,integrator,pfOrder);
  else
    return runSimulator<SIM3D>(infile,integrator,pfOrder);
}
