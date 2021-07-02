// $Id$
//==============================================================================
//!
//! \file SIMFractureDynamics.C
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for fracture-dynamic problems.
//!
//==============================================================================

#include "SIMFractureDynamics.h"

#include "SIMDynElasticity.h"
#include "SIMExplPhaseField.h"
#include "SIMPhaseField.h"
#include "SIMPoroElasticity.h"

#include "ASMbase.h"
#include "ASMunstruct.h"
#include "Functions.h"
#include "GenAlphaSIM.h"
#include "HHTSIM.h"
#include "IFEM.h"
#include "NewmarkNLSIM.h"
#include "NonLinSIM.h"
#include "SIMadmin.h"
#include "SIMCoupledSI.h"
#include "SIMenums.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TimeStep.h"
#include "Utilities.h"

#include <numeric>
#include "tinyxml.h"


bool readModel (SIMadmin& model, const std::string& infile)
{
  PROFILE("Model input");
  return model.read(infile.c_str());
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
SIMFracture<SolidSolver,PhaseSolver,Coupling>::
SIMFracture (SolidSolver& s1, PhaseSolver& s2, const std::string& inputfile)
  : CoupledSIM(s1,s2), infile(inputfile)
{
  doStop = false;
  irfStop = 0;
  E0 = Ec = Ep = stopVal = 0.0;
  refFunc = nullptr;
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
void SIMFracture<SolidSolver,PhaseSolver,Coupling>::setupDependencies ()
{
  this->S1.registerDependency(&this->S2,"phasefield",1);
  // The tensile energy is defined on integration points and not nodal points.
  // It is a global buffer array across all patches in the model.
  // Use an explicit call instead of normal couplings for this.
  this->S2.setTensileEnergy(this->S1.getTensileEnergy());
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
bool SIMFracture<SolidSolver,PhaseSolver,Coupling>::advanceStep (TimeStep& tp)
{
  if (doStop && tp.stopTime > tp.time.t)
    tp.stopTime = tp.time.t; // Update stop time due to other stop criteria
  return this->CoupledSIM::advanceStep(tp) && !doStop;
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
bool SIMFracture<SolidSolver,PhaseSolver,Coupling>::
solveStep (TimeStep& tp, bool firstS1)
{
  if (tp.step == 1 && this->S1.haveCrackPressure())
    // Start the initial step by solving the phase-field first
    if (!this->S2.solveStep(tp,false))
      return false;

  if (!this->CoupledSIM::solveStep(tp,firstS1))
    return false;

  doStop = this->S2.checkStopCriterion();
  return true;
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
bool SIMFracture<SolidSolver,PhaseSolver,Coupling>::
saveStep (const TimeStep& tp, int& nBlock)
{
  Vector RF;
  this->S1.getBoundaryReactions(RF);

  if (!energFile.empty() && tp.step > 0 &&
      this->S1.getProcessAdm().getProcId() == 0)
  {
    std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);

    Vector BF;
    this->S1.getBoundaryForce(BF,this->S1.getSolutions(),tp);

    if (tp.step == 1)
    {
      size_t i;
      os <<"#t eps_e external_energy eps+ eps- eps_b |c|"
         <<" eps_d-eps_d(0) eps_d";
      for (i = 0; i < BF.size(); i++)
        os <<" load_"<< char('X'+i);
      for (i = 0; i < RF.size(); i++)
        os <<" react_"<< char('X'+i);
      os << std::endl;
    }

    os << std::setprecision(11) << std::setw(6) << std::scientific
       << tp.time.t;
    for (double n1 : this->S1.getGlobalNorms()) os <<" "<< n1;
    const Vector& n2 = this->S2.getGlobalNorms();
    os <<" "<< (n2.size() > 2 ? n2[1] : 0.0);
    os <<" "<< (n2.size() > 1 ? n2[n2.size()-2] : 0.0);
    os <<" "<< (n2.size() > 0 ? n2.back() : 0.0);
    for (double f : BF) os <<" "<< utl::trunc(f);
    for (double f : RF) os <<" "<< utl::trunc(f);
    os << std::endl;
  }

  // Check stop criterion
  if (tp.step > 1 && irfStop > 0 && irfStop <= RF.size())
    if ((doStop = fabs(RF(irfStop)) < stopVal))
      IFEM::cout <<"\n >>> Terminating simulation due to stop criterion |RF("
                 << irfStop <<")| = "<< fabs(RF(irfStop)) <<" < "<< stopVal
                 << std::endl;

  return (this->S2.saveStep(tp,nBlock) && this->S1.saveStep(tp,nBlock) &&
          this->S2.saveResidual(tp,residual,nBlock));
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
void SIMFracture<SolidSolver,PhaseSolver,Coupling>::
parsePreref (const TiXmlElement* elem)
{
  std::string type;
  const char* value = utl::getValue(elem,elem->Value());
  if (value && utl::getAttribute(elem,"type",type))
  {
    IFEM::cout <<"\tInitial refinement function";
    refFunc = utl::parseRealFunc(value,type);
    IFEM::cout << std::endl;
  }
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
void SIMFracture<SolidSolver,PhaseSolver,Coupling>::
parseStaggering (const TiXmlElement* elem)
{
  const TiXmlElement* child = elem->FirstChildElement("stop");
  if (child)
  {
    utl::getAttribute(child,"rcomp",irfStop);
    utl::getAttribute(child,"force",stopVal);
  }
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
void SIMFracture<SolidSolver,PhaseSolver,Coupling>::
setEnergyFile (const char* fName)
{
  if (fName)
  {
    energFile = fName;
    IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
  }
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
void SIMFracture<SolidSolver,PhaseSolver,Coupling>::saveState ()
{
  sols = this->S1.getSolutions();
  sols.push_back(this->S2.getSolution());
  hsol = this->S2.getHistoryField();
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
void SIMFracture<SolidSolver,PhaseSolver,Coupling>::restoreState ()
{
  this->S1.setSolutions(Vectors(sols.begin(),sols.begin()+sols.size()-1));
  this->S2.setSolution(sols.back());
  this->S2.setHistoryField(hsol);
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
bool SIMFracture<SolidSolver,PhaseSolver,Coupling>::
preRefine (int nrefinements, int irefine, double refTol)
{
  if (irefine < 1) return true; // No pre-refinement requested

  RealFunc* refC = refFunc ? refFunc : this->S2.getInitCrack();
  if (!refC) return true; // No mesh density function defined

  // Define the initial element size on the patches
  IFEM::cout <<"\nMinimum element "
             << (SolidSolver::dimension == 3 ? "volume" : "area") <<":";
  for (const ASMbase* patch : this->S1.getFEModel())
    IFEM::cout <<" "<< patch->getMinimumSize(nrefinements);
  IFEM::cout << std::endl;

  int istat = 0, nRefine = 0;
  for (int i = 0; i < irefine && i < nrefinements; i++, refTol *= 0.5)
    if ((istat = this->S1.refine(*refC,refTol)) < 0)
      return false;
    else if (istat == 0)
      break;
    else if (!this->S2.refine(LR::RefineData()))
      return false;
    else
      ++nRefine;

  delete refFunc;
  refFunc = nullptr;
  if (nRefine == 0) return true; // No refinement

  IFEM::cout <<"\n\nPre-refinement finished. "
             <<"Reinitialize the solvers on refined mesh.\n"
             << std::string(66,'=') << std::endl;
  this->S1.clearProperties();
  this->S2.clearProperties();
  return this->S1.read(infile.c_str()) && this->S2.read(infile.c_str());
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
bool SIMFracture<SolidSolver,PhaseSolver,Coupling>::
initialRefine (double beta, double min_frac, int nrefinements)
{
  if (this->S2.hasIC("phasefield"))
    return true; // No initial refinement when specified initial phase field

  TimeStep step0;
  int newElements = 1;
  for (step0.iter = 0; newElements > 0; step0.iter++)
  {
    if (step0.iter > 0)
    {
      // Reinitialize the phase field solver (S2) on the refined mesh
      if (!readModel(this->S2,infile))
        return false;
      if (!this->S1.createFEMmodel()) // Because S2 shares the mesh with S1
        return false;
      if (!this->S2.preprocess())
        return false;
      if (!this->S2.init(step0))
        return false;
      if (!this->S2.initSystem(this->S2.opt.solver))
        return false;
    }

    IFEM::cout <<"\n\n>>> Initial refinement cycle "
               << step0.iter+1 <<" :"<< std::endl;
    if (!this->S2.solveStep(step0))
      return false;
    if ((newElements = this->adaptMesh(beta,min_frac,nrefinements,true)) < 0)
      return false;
  }

  if (step0.iter > 1)
  {
    // Reinitialize the elasticity solver (S1) on the refined mesh
    if (!readModel(this->S1,infile))
      return false;
    if (!this->S1.preprocess())
      return false;
    if (!this->S1.init(step0))
      return false;
    if (!this->S1.initSystem(this->S1.opt.solver,1,1,0,true))
      return false;
  }

  return true;
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
int SIMFracture<SolidSolver,PhaseSolver,Coupling>::
adaptMesh (double beta, double min_frac, int nrefinements,
           bool remeshOnly)
{
  // Fetch element norms to use as refinement criteria
  Vector eNorm;
  double gNorm = this->S2.getNorm(eNorm,3);
  if (eNorm.empty())
  {
    std::cerr <<" *** SIMFractureDynamics:adaptMesh: Missing refinement"
              <<" indicators, expected as the 3rd element norm."<< std::endl;
    return -1;
  }

  // Sort element indices based on comparing values in eNorm
  IntVec idx(eNorm.size());
  std::iota(idx.begin(),idx.end(),0);
  std::sort(idx.begin(),idx.end(),
            [&eNorm](size_t i1, size_t i2) { return eNorm[i1] < eNorm[i2]; });

  std::streamsize outPrec = this->S2.getOutPrec();
  std::streamsize oldPrec = outPrec > 0 ? IFEM::cout.precision(outPrec) : 0;
  double eMin = min_frac < 0.0 ? -min_frac*gNorm/sqrt(idx.size()) : min_frac;
  size_t eMax = beta < 0.0 ? idx.size() : idx.size()*beta/100.0;
  IFEM::cout <<"\n  Lowest element: "<< std::setw(8) << 1+idx.front()
             <<"    |c| = "<< eNorm[idx.front()]
             <<"\n  Highest element:"<< std::setw(8) << 1+idx.back()
             <<"    |c| = "<< eNorm[idx.back()]
             <<"\n  Minimum |c|-value for refinement: "<< eMin
             <<"\n  Minimum element "
             << (SolidSolver::dimension == 3 ? "volume" : "area") <<":";
  // Define maximum refinement level for each patch
  for (const ASMbase* patch : this->S1.getFEModel())
    IFEM::cout <<" "<< patch->getMinimumSize(nrefinements);
  IFEM::cout << std::endl;
  if (oldPrec > 0) IFEM::cout.precision(oldPrec);

  IntVec elements; // Find the elements to refine
  elements.reserve(eMax);
  for (int eid : idx)
    if (eNorm[eid] > eMin || elements.size() >= eMax)
      break;
    else
    {
      for (const ASMbase* patch : this->S1.getFEModel())
        if (patch->checkElementSize(eid))
          elements.push_back(eid);
    }

  if (elements.empty())
    return 0;

  IFEM::cout <<"  Elements to refine: "<< elements.size()
             <<" (|c| = ["<< eNorm[elements.front()]
             <<","<< eNorm[elements.back()] <<"])\n"<< std::endl;

  std::vector<LR::LRSpline*> oldBasis;
  if (!hsol.empty()) this->S2.getBasis(oldBasis);

  // Save the size of the solution vector array for the solution transfer log,
  // because refine() will resize it differently
  size_t nsol = sols.size();
  size_t nsv1 = sols.empty() ? 0 : sols.front().size();
  size_t nsv2 = sols.empty() ? 0 : sols.back().size();

  // Do the mesh refinement
  LR::RefineData prm;
  prm.options = { 10, 1, 2, 0, this->S1.getNoPatches() > 1 ? -1 : 1 };
  prm.elements = this->S1.getFunctionsForElements(elements);
  if (!this->S1.refine(prm,sols) || !this->S2.refine(prm))
    return -2;

  // Re-initialize the simulators for the new mesh
  this->S1.clearProperties();
  this->S2.clearProperties();
  if (remeshOnly)
    return elements.size();

  if (!readModel(this->S1,infile) || !readModel(this->S2,infile))
    return -3;

  if (!this->preprocess())
    return -4;

  if (!this->init(TimeStep()))
    return -5;

  if (!this->S1.initSystem(this->S1.opt.solver,1,1,0,true) ||
      !this->S2.initSystem(this->S2.opt.solver))
    return -6;

  // Transfer solution variables onto the new mesh
  if (!sols.empty())
  {
    IFEM::cout <<"\nTransferring ";
    if (nsol > 2)
      IFEM::cout << nsol-1 <<"x"<< nsv1;
    else
      IFEM::cout << nsv1;
    IFEM::cout <<" solution variables to new mesh for "<< this->S1.getName();
    Vectors soli(nsol-1,Vector(this->S1.getNoDOFs()));
    for (size_t i = 0; i < nsol-1; i++)
      for (int p = 0; p < this->S1.getNoPatches(); p++)
        this->S1.injectPatchSolution(soli[i],sols[p*nsol+i],
                                     this->S1.getPatch(p+1));
    this->S1.setSolutions(soli);
    IFEM::cout <<"\nTransferring "<< nsv2
               <<" solution variables to new mesh for "<< this->S2.getName();
    Vector solc(this->S2.getNoDOFs());
    for (int p = 0; p < this->S2.getNoPatches(); p++)
      this->S2.injectPatchSolution(solc,sols[p*nsol+nsol-1],
                                   this->S2.getPatch(p+1));
    this->S2.setSolution(solc);
  }
  if (!hsol.empty())
  {
    IFEM::cout <<"\nTransferring "<< hsol.size()
               <<" history variables to new mesh for "<< this->S2.getName()
               << std::endl;
    this->S2.transferHistory(hsol,oldBasis);
  }

  return elements.size();
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
bool SIMFracture<SolidSolver,PhaseSolver,Coupling>::
dumpMesh (const char* fileName)
{
  std::ofstream os(fileName);
  return this->S2.dumpGeometry(os);
}


template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
double SIMFracture<SolidSolver,PhaseSolver,Coupling>::
calcResidual (const TimeStep& tp, bool cycles)
{
  // Compute residual of the elasticity equation
  this->S1.setMode(SIM::RHS_ONLY);
  if (!this->S1.assembleSystem(tp.time,this->S1.getSolutions(),false))
    return -1.0;

  if (!this->S1.extractLoadVec(elastRes))
    return -1.0;

  double rNorm1 = elastRes.norm2();
  double eNorm1 = this->S1.extractScalar();

  // Compute residual of the phase-field equation
  if (!this->S2.setMode(SIM::INT_FORCES))
    return -2.0;

  Vectors sol2(1,this->S2.getSolution());
  if (!this->S2.assembleSystem(tp.time,sol2,false))
    return -2.0;

  if (!this->S2.extractLoadVec(residual))
    return -2.0;

  double rNorm2 = residual.norm2();
  double eNorm2 = this->S2.extractScalar();

  double rConv = rNorm1 + rNorm2;
  double eConv = eNorm1 + eNorm2;
  if (cycles)
  {
    IFEM::cout <<"  cycle "<< tp.iter
               <<": Res = "<< rNorm1 <<" + "<< rNorm2 <<" = "<< rConv;
    if (eConv > 0.0)
      IFEM::cout <<"  E = "<< eNorm1 <<" + "<< eNorm2 <<" = "<< eConv;
    if (tp.iter == 0)
      E0 = eConv;
    else
    {
      Ep = tp.iter > 1 ? Ec : E0;
      Ec = eConv;
      if (eConv > 0.0)
        IFEM::cout <<"  beta="<< atan2(tp.iter*(Ep-Ec),E0-Ec) * 180.0/M_PI;
    }
  }
  else
  {
    IFEM::cout <<"  Res = "<< rNorm1 <<" + "<< rNorm2 <<" = "<< rConv;
    if (eConv > 0)
      IFEM::cout <<"\n    E = "<< eNorm1 <<" + "<< eNorm2 <<" = "<< eConv;
  }
  return rConv;
}


//! \brief Helper macro to do the actual instantation.
#define INSTANCE_FULL(DIM,SIM,ELSIM,CPL) \
  template class SIMFracture<SIMDynElasticity<DIM,SIM,ELSIM<DIM>>, \
                             SIMPhaseField<DIM>,CPL>; \
  template class SIMFracture<SIMDynElasticity<DIM,SIM,ELSIM<DIM>>, \
                             SIMExplPhaseField,CPL>;

//! \brief Helper macro adding dimensionality.
#define INSTANCE_DIM(SIM,ELSIM,CPL) \
  INSTANCE_FULL(SIM2D,SIM,ELSIM,CPL) \
  INSTANCE_FULL(SIM3D,SIM,ELSIM,CPL)

//! \brief Helper macro adding nonlinear solver type.
#define INSTANCE_CPL(ELSIM,CPL) \
  INSTANCE_DIM(GenAlphaSIM,ELSIM,CPL) \
  INSTANCE_DIM(HHTSIM,ELSIM,CPL) \
  INSTANCE_DIM(NewmarkSIM,ELSIM,CPL) \
  INSTANCE_DIM(NewmarkNLSIM,ELSIM,CPL) \
  INSTANCE_DIM(NonLinSIM,ELSIM,CPL)

//! \brief Helper macro to instance for a given elasticity type.
#define INSTANCE(ELSIM) \
  INSTANCE_CPL(ELSIM,SIMCoupled) \
  INSTANCE_CPL(ELSIM,SIMCoupledSI)

INSTANCE(SIMElasticityWrap)
INSTANCE(SIMPoroElasticity)
