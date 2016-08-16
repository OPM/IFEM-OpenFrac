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

#include "InitialConditionHandler.h"
#include "CahnHilliard.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#endif
#include "TimeStep.h"
#include "Profiler.h"
#include "Utilities.h"
#include "DataExporter.h"
#include "IFEM.h"
#include "tinyxml.h"


/*!
  \brief Driver class for an Cahn-Hilliard phase-field simulator.
*/

template<class Dim> class SIMPhaseField : public Dim
{
public:
  //! \brief Default constructor.
  SIMPhaseField(Dim* gridOwner = nullptr) : Dim(1)
  {
    Dim::myHeading = "Cahn-Hilliard solver";
    if (gridOwner)
      this->clonePatches(gridOwner->getFEModel(),gridOwner->getGlob2LocMap());

    eps_d0 = refTol = 0.0;
    vtfStep = Lnorm = irefine = 0;
  }

  //! \brief Empty destructor.
  virtual ~SIMPhaseField() {}

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "CahnHilliard"; }

  //! \brief Preprocessing performed before the FEM model generation.
  virtual void preprocessA()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new CahnHilliard(Dim::dimension);

    this->printProblem();
    if (this->hasIC("phasefield"))
      IFEM::cout <<"Initial phase field specified."<< std::endl;
  }

  //! \brief Registers fields for data output.
  void registerFields(DataExporter& exporter)
  {
    exporter.registerField("c","phase field",DataExporter::SIM,
                           DataExporter::PRIMARY);
    exporter.setFieldValue("c",this,&phasefield);
  }

  //! \brief Initializes the problem.
  bool init(const TimeStep&)
  {
    this->setMode(SIM::INIT);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    this->registerField("phasefield",phasefield);
    phasefield.resize(this->getNoDOFs());
    phasefield.fill(1.0);
    return SIM::setInitialConditions(*this);
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
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMPhaseField::saveStep");

    double old = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    bool ok = this->savePoints(phasefield,tp.time.t,tp.step);
    utl::zero_print_tol = old;
    if (!ok) return false;

    if (tp.step%Dim::opt.saveInc == 0 && Dim::opt.format >= 0)
    {
      int iBlck = this->writeGlvS1(phasefield,++vtfStep,nBlock,
                                   tp.time.t,"phase",6);
      if (iBlck < 0) return false;

      // Write projected solution fields to VTF-file
      std::vector<const char*> prefix(Dim::opt.project.size(),nullptr);
      if (!Dim::opt.project.empty())
      {
        prefix.front() = Dim::opt.project.begin()->second.c_str();
        if (!this->writeGlvP(projSol,vtfStep,nBlock,iBlck,prefix.front()))
          return false;
      }

      // Write element norms to VTF-file
      const char** pfxs = prefix.empty() ? nullptr : prefix.data();
      if (!this->writeGlvN(eNorm,vtfStep,nBlock,pfxs,210))
        return false;

      if (!this->writeGlvStep(vtfStep,tp.time.t))
        return false;
    }

    if (!projSol.empty()) // Replace the phase field solution by its projection
      phasefield = projSol.getRow(1);

    return true;
  }

  //! \brief Dummy method.
  bool advanceStep(TimeStep&) { return true; }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool standalone = true)
  {
    PROFILE1("SIMPhaseField::solveStep");

    if (Dim::msgLevel == 1 && standalone)
      IFEM::cout <<"\n  Solving crack phase field at step="<< tp.step
                 <<" time="<< tp.time.t << std::endl;

    if (tp.step == 0 && tp.iter > 0) // Hack: Reduce the smearing factor
      // by a factor of 1/2 after each initial mesh refinement (at step=0)
      static_cast<CahnHilliard*>(Dim::myProblem)->scaleSmearing(0.5);

    this->setMode(SIM::STATIC);
    bool penalty = static_cast<CahnHilliard*>(Dim::myProblem)->penaltyFormulation();
    if (!this->assembleSystem(penalty ? Vectors(1,phasefield): Vectors()))
      return false;

    if (!this->solveSystem(phasefield,0))
      return false;

    if (tp.step == 1)
      static_cast<CahnHilliard*>(Dim::myProblem)->clearInitialCrack();

    return standalone ? this->postSolve(tp) : true;
  }

  //! \brief Computes solution norms, etc. on the converged solution.
  bool postSolve(TimeStep& tp)
  {
    this->printSolutionSummary(phasefield,1,
                               Dim::msgLevel > 1 ? "phasefield  " : nullptr);
    this->setMode(SIM::RECOVERY);

    // Project the phase field onto the geometry basis
    if (!Dim::opt.project.empty())
      if (!this->project(projSol,phasefield,Dim::opt.project.begin()->first))
        return false;

    Vectors gNorm;
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,Vectors(1,phasefield),gNorm,&eNorm))
      return false;
    else if (!gNorm.empty())
    {
      norm = gNorm.front();
      norm.push_back(norm.back());
      if (tp.step == 1)
        eps_d0 = norm.back();
      norm(norm.size()-1) -= eps_d0;

      if (norm.size() > 3 && utl::trunc(norm(2)) != 0.0)
      {
        if (Lnorm == 1)
          IFEM::cout <<"  L1-norm: |c^h| = (|c^h|)       : "<< norm(2)
                     <<"\n  Normalized L1-norm: |c^h|/V    : "
                     << norm(2)/norm(1) << std::endl;
        else if (Lnorm == 2)
          IFEM::cout <<"  L2-norm: |c^h| = (c^h,c^h)^0.5 : "
                     << sqrt(norm(2))
                     <<"\n  Normalized L2-norm: |c^h|/V^.5 : "
                     << sqrt(norm(2)/norm(1)) << std::endl;
      }
      if (norm.size() > 1 && utl::trunc(norm.back()) != 0.0)
        IFEM::cout <<"  Dissipated energy:       eps_d : "
                   << norm.back() << std::endl;
    }

    return true;
  }

  //! \brief Sets the tensile energy vector from the elasticity problem.
  void setTensileEnergy(const RealArray* te)
  {
    static_cast<CahnHilliard*>(Dim::myProblem)->setTensileEnergy(te);
  }

  //! \brief Returns a list of element norm values.
  double getNorm(Vector& values, size_t idx = 1) const
  {
    if (idx > eNorm.rows())
    {
      values.clear();
      return 0.0;
    }
    else if (idx == 1 || idx > 3 || Lnorm < 2)
    {
      if (Lnorm < 0 && idx > 2)
        idx--; // hack for non-present volume-specific norm
      values = eNorm.getRow(idx);
      return norm(idx);
    }

    // Apply sqrt() on the element values of the L2-norm
    values.resize(eNorm.cols());
    for (size_t i = 1; i <= values.size(); i++)
      values(i) = sqrt(eNorm(idx,i));

    return sqrt(norm(idx));
  }

  //! \brief Returns a const reference to the global norms.
  const Vector& getGlobalNorms() const { return norm; }

  //! \brief Returns a const reference to the current solution.
  const Vector& getSolution(int = 0) const { return phasefield; }

  //! \brief Updates the solution vector.
  void setSolution(const Vector& vec) { phasefield = vec; }

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
    return this->solveStep(tp,false) ? SIM::CONVERGED : SIM::DIVERGED;
  }

  //! \brief Returns the current history field.
  //! \details If projection has been done, the resulting control point values
  //! are returned, otherwise the Gauss point values are returned.
  Vector getHistoryField() const
  {
    if (projSol.rows() > 1)
      return projSol.getRow(2);

    Vector v;
    v = static_cast<const CahnHilliard*>(Dim::myProblem)->historyField;
    return v;
  }

#ifdef HAS_LRSPLINE
  //! \brief Transfers history variables at Gauss/control points to new mesh.
  //! \param[in] oldH History variables associated with Gauss- or control points
  //! \param[in] oldBasis The LR-spline basis \a oldH is referring to
  bool transferHistory2D(const RealArray& oldH, LR::LRSplineSurface* oldBasis)
  {
    const ASMu2D* pch = dynamic_cast<ASMu2D*>(this->getPatch(1));
    if (!pch) return false;

    RealArray& newH = static_cast<CahnHilliard*>(Dim::myProblem)->historyField;
    if (projSol.empty())
      return pch->transferGaussPtVars(oldBasis,oldH,newH,Dim::opt.nGauss[0]);
    else
      return pch->transferCntrlPtVars(oldBasis,oldH,newH,Dim::opt.nGauss[0]);
  }
#endif

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
        Dim::myProblem = new CahnHilliard4(Dim::dimension);
      else
        Dim::myProblem = new CahnHilliard(Dim::dimension);
    }

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"projection"))
      {
        Dim::opt.parseOutputTag(child);
        if (!Dim::opt.project.empty())
          IFEM::cout <<"\tFiltering phase field using "
                     << Dim::opt.project.begin()->second <<"."<< std::endl;
      }
      else
      {
        this->Dim::parse(child);
        // Read problem parameters (including initial crack defintition)
        if (!Dim::isRefined) // but only for the initial grid when adaptive
        {
          const char* value = utl::getValue(child,"Lnorm");
          if (value)
            Lnorm = atoi(value);
          else if ((value = utl::getValue(child,"initial_refine")))
            irefine = atoi(value);
          else if ((value = utl::getValue(child,"refine_limit")))
            refTol = atof(value);
          Dim::myProblem->parse(child);
        }
      }

#ifdef HAS_LRSPLINE
    if (Dim::isRefined)
      return true;

    // Perform initial refinement around the crack
    RealFunc* refC = static_cast<CahnHilliard*>(Dim::myProblem)->initCrack();
    ASMu2D* patch1 = dynamic_cast<ASMu2D*>(this->getPatch(1));
    if (refC && patch1)
      for (int i = 0; i < irefine; i++, refTol *= 0.5)
        if (!patch1->refine(*refC,refTol))
          return false;
#endif

    return true;
  }

private:
  Vector phasefield; //!< Current phase field solution
  Matrix projSol;    //!< Projected solution fields
  Matrix eNorm;      //!< Element norm values
  Vector norm;       //!< Global norm values
  double eps_d0;     //!< Initial eps_d value, subtracted from following values
  int    vtfStep;    //!< VTF file step counter
  int    Lnorm;      //!< Which L-norm to use to guide mesh refinement
  int    irefine;    //!< Number of initial refinement cycles
  double refTol;     //!< Initial refinement threshold
};

#endif
