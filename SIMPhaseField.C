// $Id$
//==============================================================================
//!
//! \file SIMPhaseField.C
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for Cahn-Hilliard phase-field problems.
//!
//==============================================================================

#include "SIMPhaseField.h"

#include "CahnHilliard.h"

#include "AnaSol.h"
#include "DataExporter.h"
#include "IFEM.h"
#include "Profiler.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "Vec3Oper.h"

#include "tinyxml.h"

#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "ASMu3D.h"
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#else
#include "ASMbase.h"
#endif


template<class Dim>
SIMPhaseField<Dim>::SIMPhaseField (Dim* gridOwner, size_t n) : Dim(1)
{
  Dim::myHeading = "Cahn-Hilliard solver";
  if (gridOwner && gridOwner->createFEMmodel())
  {
    this->clonePatches(gridOwner->getFEModel(),gridOwner->getGlob2LocMap());
    this->setRefined(gridOwner->getRefined());
  }
  eps_d0 = 0.0;
  vtfStep = outPrec = 0;
  transferOp = 'L';
  chp = nullptr;
  spln = nullptr;

  solution.resize(n);
}


template<class Dim>
SIMPhaseField<Dim>::~SIMPhaseField ()
{
  delete spln;
}


template<class Dim>
void SIMPhaseField<Dim>::preprocessA ()
{
  if (!Dim::myProblem)
    Dim::myProblem = chp = new CahnHilliard(Dim::dimension);

  this->printProblem();
  if (this->hasIC("phasefield"))
    IFEM::cout <<"Initial phase field specified."<< std::endl;
}


template<class Dim>
bool SIMPhaseField<Dim>::preprocessB ()
{
  if (!spln) return true;

  spln->nodes.clear();
  size_t nnod = this->getNoNodes();
  for (size_t inod = 1; inod <= nnod; inod++)
  {
    double dist = this->getNodeCoord(inod)*spln->normal + spln->d;
    if (fabs(dist) <= spln->eps) spln->nodes.push_back(inod);
  }

  IFEM::cout <<"\nSearching for nodal points on the stop plane:";
#ifdef INT_DEBUG
  for (int inod : spln->nodes)
    IFEM::cout <<"\n\t"<< inod <<"\t: "<< this->getNodeCoord(inod);
#else
  IFEM::cout <<" "<< spln->nodes.size() <<" nodes found.";
#endif
  IFEM::cout << std::endl;

  return true;
}


template<class Dim>
void SIMPhaseField<Dim>::registerFields (DataExporter& exporter)
{
  int results = DataExporter::PRIMARY;
  if (!Dim::opt.pSolOnly)
    results |= DataExporter::SECONDARY;

  exporter.registerField("c","phase field",DataExporter::SIM, results);
  exporter.setFieldValue("c",this,&solution.front());
}


template<class Dim>
bool SIMPhaseField<Dim>::init (const TimeStep& tp)
{
  bool ok = this->setMode(SIM::INIT);
  this->setQuadratureRule(Dim::opt.nGauss[0],true);
  this->registerField("phasefield",solution.front());
  solution.front().resize(this->getNoDOFs());
  if (this->hasIC("phasefield"))
    ok &= this->setInitialConditions();
  else
    solution.front().fill(1.0);
  return ok && this->advanceStep(tp);
}


template<class Dim>
bool SIMPhaseField<Dim>::saveModel (char* fileName, int& geoBlk, int& nBlock)
{
  if (Dim::opt.format < 0) return true;

  nBlock = 0;
  return this->writeGlvG(geoBlk,fileName);
}


template<class Dim>
bool SIMPhaseField<Dim>::saveStep (const TimeStep& tp, int& nBlock)
{
  PROFILE1("SIMPhaseField::saveStep");

  this->setMode(SIM::RECOVERY);

  double old = utl::zero_print_tol;
  utl::zero_print_tol = 1e-16;
  bool ok = this->savePoints(solution.front(),tp.time.t,tp.step);
  utl::zero_print_tol = old;
  if (!ok) return false;

  if (tp.step%Dim::opt.saveInc == 0 && Dim::opt.format >= 0)
  {
    int iBlck = this->writeGlvS1(solution.front(),++vtfStep,nBlock,
                                 tp.time.t,"phase",6);
    if (iBlck < 0) return false;

    // Write projected solution fields to VTF-file
    if (!Dim::opt.project.empty())
      if (!this->writeGlvP(projSol,vtfStep,nBlock,iBlck,
                           Dim::opt.project.begin()->second.c_str()))
        return false;

    if (chp->getRefinementNorm() == 2 && eNorm.rows() > 3)
      // We are using L2-norms, but the dissipated energy is not
      for (size_t i = 1; i <= eNorm.cols(); i++)
        eNorm(4,i) = eNorm(4,i)*eNorm(4,i);

    // Write element norms to VTF-file
    if (!this->writeGlvN(eNorm,vtfStep,nBlock,{},210))
      return false;

    if (!this->writeGlvStep(vtfStep,tp.time.t))
      return false;
  }

  if (transferOp == 'P') // Replace the phase field solution by its projection
    solution.front() = projSol.getRow(1);

  return true;
}


template<class Dim>
bool SIMPhaseField<Dim>::saveResidual (const TimeStep& tp,
                                       const Vector& residual, int& nBlock)
{
  if (tp.step%Dim::opt.saveInc == 0 && Dim::opt.format >= 0)
    return this->writeGlvS(residual,"phase residual",vtfStep,nBlock,110);
  else
    return true;
}


template<class Dim>
bool SIMPhaseField<Dim>::serialize (SerializeMap& data) const
{
  if (!this->saveSolution(data,this->getName()))
    return false;

  data["CH::history"] = SIMsolution::serialize(chp->historyField.data(),
                                               chp->historyField.size());
  return true;
}


template<class Dim>
bool SIMPhaseField<Dim>::deSerialize (const SerializeMap& data)
{
  if (!this->restoreSolution(data,this->getName()))
    return false;

  SerializeMap::const_iterator sit = data.find("CH::history");
  if (sit == data.end()) return false;

  SIMsolution::deSerialize(sit->second,chp->historyField.data(),
                                       chp->historyField.size());
  chp->clearInitialCrack();
  return true;
}


template<class Dim>
bool SIMPhaseField<Dim>::advanceStep (const TimeStep&)
{
  if (solution.size() > 1)
    solution.back() = solution.front();

  return true;
}


template<class Dim>
bool SIMPhaseField<Dim>::solveStep (TimeStep& tp, bool standalone)
{
  PROFILE1("SIMPhaseField::solveStep");

  if (Dim::msgLevel == 1 && standalone)
    IFEM::cout <<"\n  Solving crack phase field at step="<< tp.step
               <<" time="<< tp.time.t << std::endl;

  if (tp.step == 0 && tp.iter > 0) // Hack: Reduce the smearing factor
    // by a factor of 1/2 after each initial mesh refinement (at step=0)
    IFEM::cout <<"  - using smearing factor "<< chp->scaleSmearing(0.5)
               <<"\n"<< std::endl;

  if (!this->setMode(SIM::STATIC))
    return false;

  if (!this->assembleSystem(tp.time,solution))
    return false;

  if (!this->solveSystem(solution.front(),
                         standalone ? 0 : 1, nullptr, nullptr))
    return false;

  if (tp.step == 1)
    chp->clearInitialCrack();

  // If we solve for d as the primary phase-field variable,
  // transform it to c = 1-d here
  if (chp->useDformulation())
    for (double& c : solution.front()) c = 1.0 - c;

  return standalone ? this->postSolve(tp) : true;
}


template<class Dim>
bool SIMPhaseField<Dim>::postSolve (TimeStep& tp)
{
  this->printSolutionSummary(solution.front(), 1,
                             Dim::msgLevel > 1 ? "phasefield  " : nullptr,
                             outPrec);
  this->setMode(SIM::RECOVERY);

  // Project the phase field onto the geometry basis
  if (!Dim::opt.project.empty())
    if (!this->project(projSol,solution.front(),
                       Dim::opt.project.begin()->first))
      return false;

  Vectors gNorm;
  this->setQuadratureRule(Dim::opt.nGauss[1]);
  if (!this->solutionNorms(tp.time,solution,gNorm,&eNorm))
    return false;

  if (!gNorm.empty())
  {
    norm = gNorm.front();
    if (chp->getRefinementNorm() == 2)
      norm.back() = norm.back()*norm.back();
    norm.push_back(norm.back());
    if (tp.step == 1)
      eps_d0 = norm.back();
    norm(norm.size()-1) -= eps_d0;
    if (norm.size() > 1 && utl::trunc(norm.back()) != 0.0)
      IFEM::cout <<"  Dissipated energy:               eps_d : "
                 << norm.back() << std::endl;
    if (norm.size() > 2 && utl::trunc(norm(2)) != 0.0)
      switch (chp->getRefinementNorm())
        {
        case 1:
          IFEM::cout <<"  L1-norm: |c^h| = (|c^h|)        : "<< norm(2)
                     <<"\n  Normalized L1-norm: |c^h|/V     : "
                     << norm(2)/norm(1) << std::endl;
          break;
        case 2:
          IFEM::cout <<"  L2-norm: |c^h| = (c^h,c^h)^0.5  : "<< norm(2)
                     <<"\n  Normalized L2-norm: |c^h|/V^0.5 : "
                     << norm(2)/norm(1) << std::endl;
          break;
        }
  }
  if (gNorm.size() > 1)
  {
    static const char* norms[] = {
      "L2-norm: |c^h| = (c^h,c^h)^0.5  : ",
      "L2-norm: |c|   = (c,c)^0.5      : ",
      "L2-norm: |e|   = (e,e)^0.5      : ",
      "H1-norm: |c^h| = a(c^h,c^h)^0.5 : ",
      "H1-norm: |c|   = a(c,c)^0.5     : ",
      "H1-norm: |e|   = a(e,e)^0.5     : ",
      nullptr };
    const Vector& nrm = gNorm.back();
    for (size_t i = 0; i < nrm.size() && norms[i]; i++)
      if (utl::trunc(nrm[i]) != 0.0)
        IFEM::cout <<"\n  "<< norms[i] << nrm[i];
    IFEM::cout << std::endl;
  }

#ifdef HAVE_MPI
  if (this->adm.dd.isPartitioned())
    this->adm.allReduce(chp->historyField, MPI_MAX);
#endif

  return true;
}


template<class Dim>
void SIMPhaseField<Dim>::printSolutionSummary (const Vector& solvec,
                                               int printSol,
                                               const char* compName,
                                               std::streamsize prec)
{
  this->Dim::printSolutionSummary(solvec,printSol,compName,prec);

  // Also print the smallest phase field value and the range
  int inod = 0, minNod = 0;
  double minVal = 1.0, maxVal = 1.0;
  for (double v : solution.front())
  {
    ++inod;
    if (v < minVal)
    {
      minVal = v;
      minNod = inod;
    }
    else if (v > maxVal)
      maxVal = v;
  }

  std::streamsize oldPrec = prec > 0 ? IFEM::cout.precision(prec) : 0;
  IFEM::cout <<"                            Min value    : "
             << minVal <<" node "<< minNod
             <<"\n                            Range        : "
             << maxVal-minVal << std::endl;
  if (oldPrec > 0) IFEM::cout.precision(oldPrec);
}


template<class Dim>
bool SIMPhaseField<Dim>::checkStopCriterion () const
{
  if (!spln || spln->nodes.empty())
    return false;

  int imin = 0;
  double c = 0.0, cmin = 1.0;
  for (int inod : spln->nodes)
    if ((c = solution.front()[inod-1]) < spln->cstop)
    {
      IFEM::cout <<"\n >>> Terminating simulation due to stop criterion c("
                 << this->getNodeCoord(inod) <<") = "<< c <<" < "
                 << spln->cstop << std::endl;
      return true;
    }
    else if (c < cmin)
    {
      cmin = c;
      imin = inod;
    }

  IFEM::cout <<"  Minimum phase field value on plane: "<< cmin <<" at node "
             << imin <<" (X = "<< this->getNodeCoord(imin) <<")"<< std::endl;
  return false;
}


template<class Dim>
RealFunc* SIMPhaseField<Dim>::getInitCrack() const
{
  return chp->initCrack();
}


template<class Dim>
void SIMPhaseField<Dim>::setTensileEnergy (const RealArray* te)
{
  chp->setTensileEnergy(te);
}


template<class Dim>
double SIMPhaseField<Dim>::getNorm (Vector& values, size_t idx) const
{
  if (idx < 1 || idx > eNorm.rows())
  {
    values.clear();
    return 0.0;
  }
  else
  {
    values = eNorm.getRow(idx);
    return norm(idx);
  }
}


template<class Dim>
SIM::ConvStatus SIMPhaseField<Dim>::solveIteration (TimeStep& tp)
{
  if (Dim::msgLevel == 1)
    IFEM::cout <<"\n  step="<< tp.step <<"  time="<< tp.time.t << std::endl;
  return this->solveStep(tp,false) ? SIM::CONVERGED : SIM::FAILURE;
}


template<class Dim>
RealArray SIMPhaseField<Dim>::getHistoryField () const
{
  if (transferOp == 'P' && projSol.rows() > 1)
    return projSol.getRow(2);

  return chp->historyField;
}


template<class Dim>
void SIMPhaseField<Dim>::setHistoryField (const RealArray& hfield)
{
  RealArray& hist = chp->historyField;
  if (transferOp != 'P' && hist.size() == hfield.size())
    std::copy(hfield.begin(),hfield.end(),hist.begin());
}


template<class Dim>
void SIMPhaseField<Dim>::getBasis (std::vector<LR::LRSpline*>& basis)
{
#ifdef HAS_LRSPLINE
  for (ASMbase* patch : Dim::myModel)
  {
    ASMu2D* pch2 = dynamic_cast<ASMu2D*>(patch);
    if (pch2)
      basis.push_back(pch2->getBasis()->copy());
    else
    {
      ASMu3D* pch3 = dynamic_cast<ASMu3D*>(patch);
      if (pch3)
        basis.push_back(pch3->getBasis()->copy());
      else
        continue;
    }
    basis.back()->generateIDs();
  }
#endif
}


template<class Dim>
bool SIMPhaseField<Dim>::transferHistory (const RealArray& oldH,
                                          std::vector<LR::LRSpline*>& oldBasis)
{
  RealArray& newH = chp->historyField;
  newH.clear();

#ifdef HAS_LRSPLINE
  bool ok = true;
  int nGp = Dim::opt.nGauss[0];
  RealArray::const_iterator itH = oldH.begin(), jtH;
  RealArray newHp, oldHn;
  size_t idx = 0;
  for (LR::LRSpline* basis : oldBasis)
  {
    ASMLRSpline* pch = dynamic_cast<ASMLRSpline*>(this->getPatch(++idx));
    if (pch && idx <= oldBasis.size())
    {
      switch (transferOp) {
      case 'P':
        oldHn.resize(pch->getNoNodes());
        for (size_t i = 0; i < oldHn.size(); i++)
          oldHn[i] = oldH[pch->getNodeID(1+i)-1];
        ok &= pch->transferCntrlPtVars(basis,oldHn,newHp,nGp);
        break;
      case 'N':
        jtH = itH + basis->nElements()*pow(nGp, static_cast<int>(Dim::dimension));
        ok &= pch->transferGaussPtVarsN(basis,RealArray(itH,jtH),newHp,nGp);
        itH = jtH;
        break;
      default:
        jtH = itH + basis->nElements()*pow(nGp, static_cast<int>(Dim::dimension));
        ok &= pch->transferGaussPtVars(basis,RealArray(itH,jtH),newHp,nGp);
        itH = jtH;
      }
      newH.insert(newH.end(),newHp.begin(),newHp.end());
    }
    delete basis;
  }

  return ok;
#else
  return false;
#endif
}


template<class Dim>
bool SIMPhaseField<Dim>::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"postprocessing"))
    utl::getAttribute(elem,"precision",outPrec);

  if (strcasecmp(elem->Value(),"cahnhilliard"))
    return this->Dim::parse(elem);

  if (!Dim::myProblem)
  {
    int order = 2;
    utl::getAttribute(elem,"order",order);
    if (order == 4)
      Dim::myProblem = chp = new CahnHilliard4(Dim::dimension);
    else
      Dim::myProblem = chp = new CahnHilliard(Dim::dimension);
  }

  bool result = true;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"projection"))
    {
      Dim::opt.parseOutputTag(child);
      if (!Dim::opt.project.empty())
        IFEM::cout <<"\tFiltering phase field using "
                   << Dim::opt.project.begin()->second <<"."<< std::endl;
    }
    else if (!strncasecmp(child->Value(),"use_project",11) &&
             !Dim::opt.project.empty())
    {
      transferOp = 'P';
      IFEM::cout <<"\tUsing projected solution transfer."<< std::endl;
    }
    else if (!strcasecmp(child->Value(),"use_NNtransfer"))
    {
      transferOp = 'N';
      IFEM::cout <<"\tUsing nearest-neighbour solution transfer."<< std::endl;
    }
    else if (!strcasecmp(child->Value(),"anasol"))
    {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type.empty())
        type = "Expression";
      else
        type[0] = toupper(type[0]);
      IFEM::cout <<"\tAnalytical solution: "<< type << std::endl;
      if (!Dim::mySol)
        Dim::mySol = new AnaSol(child,true);
    }
    else if (!strcasecmp(child->Value(),"stop_plane"))
    {
      if (!spln) spln = new StopPlane();
      utl::getAttribute(child,"nx",spln->normal.x);
      utl::getAttribute(child,"ny",spln->normal.y);
      utl::getAttribute(child,"nz",spln->normal.z);
      utl::getAttribute(child,"d",spln->d);
      utl::getAttribute(child,"eps",spln->eps);
      utl::getAttribute(child,"stop",spln->cstop);
      IFEM::cout <<"\tStop plane: normal="<< spln->normal
                 <<" d="<< spln->d
                 <<"\n\t            eps="<< spln->eps <<" stop="<< spln->cstop
                 << std::endl;
    }
    else
    {
      if (!strcasecmp(child->Value(),"postprocessing"))
        utl::getAttribute(child,"precision",outPrec);
      result &= this->Dim::parse(child);
      // Read problem parameters (including initial crack definition)
      chp->parse(child,Dim::isRefined);
    }

  return result;
}


template<class Dim>
bool SIMPhaseField<Dim>::initNeumann (size_t propInd)
{
  typename Dim::VecFuncMap::const_iterator tit = Dim::myVectors.find(propInd);
  if (tit == Dim::myVectors.end()) return false;

  chp->setFlux(tit->second);

  return true;
}


template class SIMPhaseField<SIM1D>;
template class SIMPhaseField<SIM2D>;
template class SIMPhaseField<SIM3D>;
