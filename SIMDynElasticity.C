// $Id$
//==============================================================================
//!
//! \file SIMDynElasticity.C
//!
//! \date Dec 04 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dynamic simulation driver for elasticity problems with fracture.
//!
//==============================================================================

#include "SIMDynElasticity.h"

#include "FractureElasticityVoigt.h"
#include "SIMPoroElasticity.h"

#include "GenAlphaSIM.h"
#include "HHTSIM.h"
#include "NewmarkSIM.h"
#include "NewmarkNLSIM.h"
#include "NonLinSIM.h"
#include "SIM2D.h"
#include "SIM3D.h"

#ifdef IFEM_HAS_POROELASTIC
#include "PoroFracture.h"
#endif

#include <fstream>


template< class Dim, class DynSIM, class Sim>
SIMDynElasticity<Dim,DynSIM,Sim>::SIMDynElasticity () :
  dSim(*this)
{
  vtfStep = outPrec = 0;
}


template< class Dim, class DynSIM, class Sim>
SIMDynElasticity<Dim,DynSIM,Sim>::
SIMDynElasticity(const std::vector<unsigned char>& nf)
  : Sim(nf), dSim(*this)
{
  vtfStep = outPrec = 0;
}


template< class Dim, class DynSIM, class Sim>
void SIMDynElasticity<Dim,DynSIM,Sim>::printProblem () const
{
  static short int ncall = 0;
  if (++ncall == 1) // Avoiding infinite recursive calls
    dSim.printProblem();
  else
    this->Sim::printProblem();
  --ncall;
}


template< class Dim, class DynSIM, class Sim>
bool SIMDynElasticity<Dim,DynSIM,Sim>::init (const TimeStep& tp, bool)
{
  dSim.initPrm();
  dSim.initSol(dynamic_cast<NewmarkSIM*>(&dSim) ? 3 : 2);

  bool ok = this->setMode(SIM::INIT) && this->getIntegrand()->init(tp.time);
  this->setQuadratureRule(Dim::opt.nGauss[0],true);
  this->registerField("solution",dSim.getSolution());
  return this->setInitialConditions() && ok;
}


template< class Dim, class DynSIM, class Sim>
bool SIMDynElasticity<Dim,DynSIM,Sim>::saveStep (const TimeStep& tp,
                                                 int& nBlock)
{
  double old = utl::zero_print_tol;
  utl::zero_print_tol = 1e-16;
  bool ok = this->savePoints(dSim.getSolution(),tp.time.t,tp.step);
  utl::zero_print_tol = old;

  if (!energFile.empty() && tp.step > 0 && Dim::adm.getProcId() == 0)
  {
    std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);

    Vector Bforce, Rforce;
    this->getBoundaryForce(Bforce,dSim.getSolutions(),tp);
    this->getBoundaryReactions(Rforce);

    if (tp.step == 1)
    {
      size_t i;
      os <<"#t eps_e external_energy eps+ eps- eps_b";
      for (i = 0; i < Bforce.size(); i++)
        os <<" load_"<< char('X'+i);
      for (i = 0; i < Rforce.size(); i++)
        os <<" react_"<< char('X'+i);
      os << std::endl;
    }

    os << std::setprecision(11) << std::setw(6) << std::scientific
       << tp.time.t;
    for (double n : gNorm)  os <<" "<< n;
    for (double f : Bforce) os <<" "<< utl::trunc(f);
    for (double f : Rforce) os <<" "<< utl::trunc(f);
    os << std::endl;
  }

  if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0 || !ok)
    return ok;

  // Write primary and secondary (of requested) solution fields to VTF-file
  if (!dSim.saveStep(++vtfStep,nBlock,tp.time.t))
    return false;
  else if (tp.step < 1)
    return true;

  // Write projected solution fields to VTF-file
  if (!Dim::opt.project.empty())
    if (!this->writeGlvP(projSol,vtfStep,nBlock,110,
                         Dim::opt.project.begin()->second.c_str()))
      return false;

  // Write element norms to VTF-file
  return this->writeGlvN(eNorm,vtfStep,nBlock);
}


template< class Dim, class DynSIM, class Sim>
void SIMDynElasticity<Dim,DynSIM,Sim>::printStep (int istep,
                                                  const TimeDomain& time) const
{
  Dim::adm.cout <<"\n  step="<< istep <<"  time="<< time.t;
  RealFunc* p = this->haveCrackPressure();
  if (p && !p->isConstant())
    Dim::adm.cout <<"  crack pressure: "
                  << (*p)(Vec4(0.0,0.0,0.0,time.t));
  Dim::adm.cout << std::endl;
}


template< class Dim, class DynSIM, class Sim>
bool SIMDynElasticity<Dim,DynSIM,Sim>::solveStep (TimeStep& tp)
{
  if (Dim::msgLevel >= 1)
    IFEM::cout <<"\n  Solving the elasto-dynamics problem...";

  if (dSim.solveStep(tp,SIM::STATIC,1.0e-8,outPrec) != SIM::CONVERGED)
    return false;

  return this->postSolve(tp);
}


template< class Dim, class DynSIM, class Sim>
bool SIMDynElasticity<Dim,DynSIM,Sim>::postSolve (TimeStep& tp)
{
  RealArray RF;
  if (this->getCurrentReactions(RF,dSim.getSolution()))
  {
    IFEM::cout <<"  Total reaction forces:          Sum(R) :";
    for (size_t i = 1; i < RF.size(); i++)
      IFEM::cout <<" "<< utl::trunc(RF[i]);
    double Ru = RF.front();
    if (utl::trunc(Ru) != 0.0)
      IFEM::cout <<"\n  displacement*reactions:          (R,u) : "<< Ru;
    IFEM::cout << std::endl;
  }

#ifdef HAVE_MPI
  if (this->adm.dd.isPartitioned()) {
    RealArray& hist = const_cast<RealArray&>(*this->getTensileEnergy());
    std::fill(hist.begin(), hist.end(), 0.0);
  }
#endif

  // Update strain energy density for the converged solution
  if (!this->updateStrainEnergyDensity(tp))
    return false;

#ifdef HAVE_MPI
  if (this->adm.dd.isPartitioned()) {
    RealArray& hist = const_cast<RealArray&>(*this->getTensileEnergy());
    this->adm.allReduceAsSum(hist);
  }
#endif
  // Project the secondary solution field onto the geometry basis
  if (!Dim::opt.project.empty())
    if (!this->project(projSol,dSim.getSolution(),
                       Dim::opt.project.begin()->first))
      return false;

  Vectors gNorms;
  this->setQuadratureRule(Dim::opt.nGauss[1]);
  if (!this->solutionNorms(tp.time,dSim.getSolutions(),gNorms,&eNorm))
    return false;
  else if (!gNorms.empty())
  {
    std::streamsize oldPrec = outPrec > 0 ? IFEM::cout.precision(outPrec) : 0;
    gNorm = gNorms.front();
    if (gNorm.size() > 0 && utl::trunc(gNorm(1)) != 0.0)
      IFEM::cout <<"  Elastic strain energy:           eps_e : "
                 << gNorm(1) << std::endl;
    if (gNorm.size() > 4 && utl::trunc(gNorm(5)) != 0.0)
      IFEM::cout <<"  Bulk energy:                     eps_b : "
                 << gNorm(5)
                 <<"\n  Tensile & compressive energies         : "
                 << gNorm(3) <<" "<< gNorm(4) << std::endl;
    if (gNorm.size() > 1 && utl::trunc(gNorm(2)) != 0.0)
      IFEM::cout <<"  External energy: ((f,u^h)+(t,u^h))^0.5 : "
                 << (gNorm(2) < 0.0 ? -sqrt(-gNorm(2)) : sqrt(gNorm(2)))
                 << std::endl;
    if (oldPrec > 0) IFEM::cout.precision(oldPrec);
  }

  if (this->hasResultPoints())
    this->dumpResults(dSim.getSolution(), tp.time.t, IFEM::cout, true, 5);

  return true;
}


template< class Dim, class DynSIM, class Sim>
bool SIMDynElasticity<Dim,DynSIM,Sim>::
updateStrainEnergyDensity (const TimeStep& tp)
{
  this->setMode(SIM::RECOVERY);
  return this->assembleSystem(tp.time,dSim.getSolutions());
}


template< class Dim, class DynSIM, class Sim>
const RealArray* SIMDynElasticity<Dim,DynSIM,Sim>::getTensileEnergy () const
{
  return static_cast<Elasticity*>(Dim::myProblem)->getTensileEnergy();
}


template< class Dim, class DynSIM, class Sim>
void SIMDynElasticity<Dim,DynSIM,Sim>::setEnergyFile (const char* fName)
{
  if (fName)
  {
    energFile = fName;
    IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
  }
}


template< class Dim, class DynSIM, class Sim>
void SIMDynElasticity<Dim,DynSIM,Sim>::setSolutions (const Vectors& dvec)
{
  size_t nSol = dSim.getSolutions().size();
  for (size_t i = 0; i < nSol && i < dvec.size(); i++)
    dSim.setSolution(dvec[i],i);
}


template< class Dim, class DynSIM, class Sim>
SIM::ConvStatus SIMDynElasticity<Dim,DynSIM,Sim>::
solveIteration(TimeStep& tp, char stage)
{
  if (Dim::msgLevel == 1 && tp.iter == 0 && stage == 1)
    this->printStep(tp.step,tp.time);
  dSim.setSubIteration(tp.iter == 0 ? DynSIM::FIRST : DynSIM::ITER);

  if (tp.iter == 0 && stage == 1)
    return dSim.solveIteration(tp);
  else if (tp.iter > 0 && stage == 2)
  {
    SIM::ConvStatus status = SIM::OK;
    while (tp.iter <= dSim.getMaxit())
      switch ((status = dSim.solveIteration(tp))) {
      case SIM::OK:
      case SIM::SLOW:
        tp.iter++;
        break;
      default:
        return status;
      }
    return SIM::DIVERGED; // No convergence in maxit iterations
  }

  // Solve the whole time step
  TimeStep myTp(tp); // Make a copy to avoid destroying the iteration counter
  return dSim.solveStep(myTp);
}



template< class Dim, class DynSIM, class Sim>
RealFunc* SIMDynElasticity<Dim,DynSIM,Sim>::haveCrackPressure () const
{
#ifdef IFEM_HAS_POROELASTIC
  PoroFracture* pfel = dynamic_cast<PoroFracture*>(Dim::myProblem);
  if (pfel) return pfel->getCrackPressure();
#endif
  FractureElasticity* fel = dynamic_cast<FractureElasticity*>(Dim::myProblem);
  return fel ? fel->getCrackPressure() : nullptr;
}


template< class Dim, class DynSIM, class Sim>
bool SIMDynElasticity<Dim,DynSIM,Sim>::
serialize (SIMsolution::SerializeMap& data) const
{
  if (!this->SIMElasticityWrap<Dim>::serialize(data))
    return false;

  data["DynSIM::refNorm"] = SIMsolution::serialize(dSim.getRefNorm(),1);

  return true;
}


template< class Dim, class DynSIM, class Sim>
bool SIMDynElasticity<Dim,DynSIM,Sim>::
deSerialize(const SIMsolution::SerializeMap& data)
{
  if (!this->SIMElasticityWrap<Dim>::deSerialize(data))
    return false;

  SIMsolution::SerializeMap::const_iterator it;
  if ((it = data.find("DynSIM::refNorm")) != data.end())
    SIMsolution::deSerialize(it->second,dSim.theRefNorm(),1);

  return true;
}


template< class Dim, class DynSIM, class Sim>
Elasticity* SIMDynElasticity<Dim,DynSIM,Sim>::getIntegrand ()
{
  if (!Dim::myProblem) // Using the Voigt formulation by default
    Dim::myProblem = new FractureElasticityVoigt(Dim::dimension);
  return static_cast<Elasticity*>(Dim::myProblem);
}


template< class Dim, class DynSIM, class Sim>
bool SIMDynElasticity<Dim,DynSIM,Sim>::parse (const TiXmlElement* elem)
{
  bool result = true;
  static short int ncall = 0;
  if (++ncall == 1) // Avoiding infinite recursive calls
    result = dSim.parse(elem);
  else if (!strcasecmp(elem->Value(),SIMElasticity<Dim>::myContext.c_str()))
  {
    if (!Dim::myProblem)
    {
      if (this->getName() == "PoroElasticity")
#ifdef IFEM_HAS_POROELASTIC
        Dim::myProblem = new PoroFracture(Dim::dimension,
                                          this->mixedProblem());
#else
        return false;
#endif
      else
      {
        std::string formulation("voigt");
        utl::getAttribute(elem,"formulation",formulation,true);
        if (formulation != "voigt")
          Dim::myProblem = new FractureElasticity(Dim::dimension);
      }
    }
    result = this->Sim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"noGeometricStiffness"))
        this->setIntegrationPrm(3,1); // Disable geometric stiffness
      else if (!strcasecmp(child->Value(),"postprocessing"))
        utl::getAttribute(child,"precision",outPrec);
  }
  else
  {
    if (!strcasecmp(elem->Value(),"postprocessing"))
      utl::getAttribute(elem,"precision",outPrec);
    result = this->Dim::parse(elem);
  }
  --ncall;
  return result;
}


//! \brief Helper macro to perform actual instantation.
#define INSTANCE_DIM(SIM,ELSIM) \
  template class SIMDynElasticity<SIM2D,SIM,ELSIM<SIM2D>>; \
  template class SIMDynElasticity<SIM3D,SIM,ELSIM<SIM3D>>;


//! \brief Helper macro instanting for elasticity variants.
#define INSTANCE(SIM) \
  INSTANCE_DIM(SIM,SIMPoroElasticity) \
  INSTANCE_DIM(SIM,SIMElasticityWrap)

INSTANCE(GenAlphaSIM)
INSTANCE(HHTSIM)
INSTANCE(LinSIM)
INSTANCE(NewmarkSIM)
INSTANCE(NewmarkNLSIM)
INSTANCE(NonLinSIM)
