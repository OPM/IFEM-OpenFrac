// $Id$
//==============================================================================
//!
//! \file SIMPhaseField.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for Cahn-Hilliard phase-field problems.
//!
//==============================================================================

#ifndef _SIM_PHASE_FIELD_H
#define _SIM_PHASE_FIELD_H

#include "CahnHilliard.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "ASMu3D.h"
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#else
namespace LR { class LRSpline; }
#endif
#include "TimeStep.h"
#include "Profiler.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "DataExporter.h"
#include "IFEM.h"
#include "tinyxml.h"


/*!
  \brief Driver class for a Cahn-Hilliard phase-field simulator.
*/

template<class Dim> class SIMPhaseField : public Dim
{
  //! \brief Struct defining a user-specified plane for stopping the simulation.
  struct StopPlane
  {
    Vec3             normal; //!< Normal vector of the stop plane
    double           d;      //!< Distance offset of plane origin
    double           eps;    //!< Zero tolerance for nodes being on the plane
    std::vector<int> nodes;  //!< Nodal (control) points on the plane
    double           cstop;  //!< Stop when phase field value lower that this
  };

public:
  //! \brief Default constructor.
  explicit SIMPhaseField(Dim* gridOwner = nullptr, size_t n = 2) : Dim(1)
  {
    Dim::myHeading = "Cahn-Hilliard solver";
    if (gridOwner && gridOwner->createFEMmodel())
    {
      Dim::nf.resize(gridOwner->getNoBasis(),0);
      this->clonePatches(gridOwner->getFEModel(),gridOwner->getGlob2LocMap());
    }

    eps_d0 = refTol = 0.0;
    vtfStep = irefine = 0;
    transferOp = 'L';
    chp = nullptr;
    spln = nullptr;

    phasefield.resize(n);
  }

  //! \brief The destructor deletes the stop plane.
  virtual ~SIMPhaseField() { delete spln; }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "CahnHilliard"; }

  //! \brief Preprocessing performed before the FEM model generation.
  virtual void preprocessA()
  {
    if (!Dim::myProblem)
      Dim::myProblem = chp = new CahnHilliard(Dim::dimension);

    this->printProblem();
    if (this->hasIC("phasefield"))
      IFEM::cout <<"Initial phase field specified."<< std::endl;
  }

  //! \brief Preprocessing performed after the FEM model generation.
  virtual bool preprocessB()
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

  //! \brief Registers fields for data output.
  void registerFields(DataExporter& exporter)
  {
    exporter.registerField("c","phase field",DataExporter::SIM,
                           DataExporter::PRIMARY);
    exporter.setFieldValue("c",this,&phasefield.front());
  }

  //! \brief Initializes the problem.
  bool init(const TimeStep& tp)
  {
    bool ok = this->setMode(SIM::INIT);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    this->registerField("phasefield",phasefield.front());
    if (this->hasIC("phasefield"))
      for (Vector& solvec : phasefield)
        solvec.resize(this->getNoDOFs(),1.0);
    return this->setInitialConditions() && ok && this->advanceStep(tp);
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMPhaseField::saveStep");

    double old = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    bool ok = this->savePoints(phasefield.front(),tp.time.t,tp.step);
    utl::zero_print_tol = old;
    if (!ok) return false;

    if (tp.step%Dim::opt.saveInc == 0 && Dim::opt.format >= 0)
    {
      int iBlck = this->writeGlvS1(phasefield.front(),++vtfStep,nBlock,
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
      if (!this->writeGlvN(eNorm,vtfStep,nBlock,nullptr,210))
        return false;

      if (!this->writeGlvStep(vtfStep,tp.time.t))
        return false;
    }

    if (transferOp == 'P') // Replace the phase field solution by its projection
      phasefield.front() = projSol.getRow(1);

    return true;
  }

  //! \brief Saves the force residual of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] residual Residual force vector
  //! \param nBlock Running VTF block counter
  bool saveResidual(const TimeStep& tp, const Vector& residual, int& nBlock)
  {
    if (tp.step%Dim::opt.saveInc == 0 && Dim::opt.format >= 0)
      return this->writeGlvS(residual,"phase residual",vtfStep,nBlock,110);
    else
      return true;
  }

  //! \brief Advances the time step one step forward.
  bool advanceStep(const TimeStep&)
  {
    if (phasefield.size() > 1)
      phasefield.back() = phasefield.front();
    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool standalone = true)
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

    if (!this->assembleSystem(tp.time,phasefield))
      return false;

    if (!this->solveSystem(phasefield.front(),
                           standalone ? 0 : 1, nullptr, nullptr))
      return false;

    if (tp.step == 1)
      chp->clearInitialCrack();

    // If we solve for d as the primary phase-field variable,
    // transform it to c = 1-d here
    if (chp->useDformulation())
      for (double& c : phasefield.front()) c = 1.0 - c;

    return standalone ? this->postSolve(tp) : true;
  }

  //! \brief Computes solution norms, etc. on the converged solution.
  bool postSolve(TimeStep& tp)
  {
    this->printSolutionSummary(phasefield.front(),1,
                               Dim::msgLevel > 1 ? "phasefield  " : nullptr);
    this->setMode(SIM::RECOVERY);

    // Project the phase field onto the geometry basis
    if (!Dim::opt.project.empty())
      if (!this->project(projSol,phasefield.front(),
                         Dim::opt.project.begin()->first))
        return false;

    Vectors gNorm;
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,phasefield,gNorm,&eNorm))
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

    return true;
  }

  //! \brief Prints a summary of the calculated solution to console.
  //! \param[in] solution The solution vector
  //! \param[in] printSol Print solution only if size is less than this value
  //! \param[in] compName Solution name to be used in norm output
  virtual void printSolutionSummary (const Vector& solution, int printSol,
                                     const char* compName, std::streamsize = 0)
  {
    this->Dim::printSolutionSummary(solution,printSol,compName);

    // Also print the smallest phase field value and the range
    int inod = 0, minNod = 0;
    double minVal = 1.0, maxVal = 1.0;
    for (double v : phasefield.front())
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

    if (minNod > 0)
      IFEM::cout <<"                            Min value    : "
                 << minVal <<" node "<< minNod
                 <<"\n                            Range        : "
                 << maxVal-minVal << std::endl;
  }

  //! \brief Returns \e true if terminating due to user-defined criteria.
  bool checkStopCriterion () const
  {
    if (!spln || spln->nodes.empty())
      return false;

    int imin = 0;
    double c = 0.0, cmin = 1.0;
    for (int inod : spln->nodes)
      if ((c = phasefield.front()[inod-1]) < spln->cstop)
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

  //! \brief Sets the tensile energy vector from the elasticity problem.
  void setTensileEnergy(const RealArray* te) { chp->setTensileEnergy(te); }

  //! \brief Returns a list of element norm values.
  double getNorm(Vector& values, size_t idx = 1) const
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

  //! \brief Returns a const reference to the global norms.
  const Vector& getGlobalNorms() const { return norm; }

  //! \brief Returns a const reference to the current solution.
  const Vector& getSolution() const { return phasefield.front(); }
  //! \brief Updates the solution vector.
  void setSolution(const Vector& vec) { phasefield.front() = vec; }

  //! \brief Returns the maximum number of iterations (unlimited).
  int getMaxit() const { return 9999; }
  //! \brief Returns the number of initial refinement cycles.
  int getInitRefine() const { return irefine; }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //!
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    if (Dim::msgLevel == 1)
      IFEM::cout <<"\n  step="<< tp.step <<"  time="<< tp.time.t << std::endl;
    return this->solveStep(tp,false) ? SIM::CONVERGED : SIM::FAILURE;
  }

  //! \brief Returns the current history field.
  //! \details If projection has been done, the resulting control point values
  //! are returned, otherwise the Gauss point values are returned.
  RealArray getHistoryField() const
  {
    if (transferOp == 'P' && projSol.rows() > 1)
      return projSol.getRow(2);

    return chp->historyField;
  }

  //! \brief Resets the history field from the provided array.
  void setHistoryField(const RealArray& hfield)
  {
    RealArray& hist = chp->historyField;
    if (transferOp != 'P' && hist.size() == hfield.size())
      std::copy(hfield.begin(),hfield.end(),hist.begin());
  }

  //! \brief Extracts the LR-spline basis for the phase field
  void getBasis(std::vector<LR::LRSpline*>& basis)
  {
#ifdef HAS_LRSPLINE
    ASMu2D* pch2;
    ASMu3D* pch3;
    for (ASMbase* patch : Dim::myModel)
      if ((pch2 = dynamic_cast<ASMu2D*>(patch)))
        basis.push_back(pch2->getBasis()->copy());
      else if ((pch3 = dynamic_cast<ASMu3D*>(patch)))
        basis.push_back(pch3->getBasis()->copy());
#endif
  }

  //! \brief Transfers history variables at Gauss/control points to new mesh.
  //! \param[in] oldH History variables associated with Gauss- or control points
  //! \param[in] oldBasis The LR-spline basis \a oldH is referring to
  bool transferHistory(const RealArray& oldH,
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
      ASMunstruct* pch = dynamic_cast<ASMunstruct*>(this->getPatch(++idx));
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

protected:
  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
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
        result &= this->Dim::parse(child);
        // Read problem parameters (including initial crack defintition)
        if (!Dim::isRefined) // but only for the initial grid when adaptive
        {
          const char* value = utl::getValue(child,"initial_refine");
          if (value)
            irefine = atoi(value);
          else if ((value = utl::getValue(child,"refine_limit")))
            refTol = atof(value);
          Dim::myProblem->parse(child);
        }
      }

#ifdef HAS_LRSPLINE
    if (Dim::isRefined || !result)
      return result;

    // Perform initial refinement around the crack (single-patch only)
    RealFunc* refC = chp->initCrack();
    ASMunstruct* patch1 = dynamic_cast<ASMunstruct*>(this->getPatch(1));
    if (refC && patch1 && this->getNoPatches() == 1)
      for (int i = 0; i < irefine; i++, refTol *= 0.5)
        if (!patch1->refine(*refC,refTol))
          return false;
#endif

    return result;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::VecFuncMap::const_iterator tit = Dim::myVectors.find(propInd);
    if (tit == Dim::myVectors.end()) return false;

    chp->setFlux(tit->second);

    return true;
  }

private:
  CahnHilliard* chp;  //!< The Cahn-Hilliard integrand
  StopPlane*    spln; //!< Stop simulation when crack penetrates this plane

  Vectors phasefield; //!< Current (and previous) phase field solution
  Matrix  projSol;    //!< Projected solution fields
  Matrix  eNorm;      //!< Element norm values
  Vector  norm;       //!< Global norm values
  double  eps_d0;     //!< Initial eps_d value, subtracted from following values
  int     vtfStep;    //!< VTF file step counter
  int     irefine;    //!< Number of initial refinement cycles
  double  refTol;     //!< Initial refinement threshold
  char    transferOp; //!< Solution transfer option
};

#endif
