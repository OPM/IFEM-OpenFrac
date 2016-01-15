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
#include "TimeStep.h"
#include "Profiler.h"
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
  SIMPhaseField(char order = 2) : Dim(1), eps_d0(0.0)
  {
    Dim::myHeading = "Cahn-Hilliard solver";
    if (order == 4)
      Dim::myProblem = new CahnHilliard4(Dim::dimension);
    else
      Dim::myProblem = new CahnHilliard(Dim::dimension);
  }

  //! \brief Empty destructor.
  virtual ~SIMPhaseField() {}

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "CahnHilliard"; }

  //! \brief Preprocessing performed before the FEM model generation.
  virtual void preprocessA()
  {
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
    return true;
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

    if (tp.step%Dim::opt.saveInc == 0 && Dim::opt.format >= 0)
    {
      int iBlck = 6;
      int iDump = tp.step/Dim::opt.saveInc;
      iBlck = this->writeGlvS1(phasefield,iDump,nBlock,tp.time.t,"phase",iBlck);
      if (iBlck < 0) return false;

      if (!Dim::opt.project.empty())
        if (!this->writeGlvP(projSol,iDump,nBlock,iBlck,
                             Dim::opt.project.begin()->second.c_str()))
          return false;

      if (!this->writeGlvStep(iDump,tp.time.t))
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

    this->setMode(SIM::STATIC);
    if (!this->assembleSystem())
      return false;

    if (!this->solveSystem(phasefield,0))
      return false;

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
    if (!this->solutionNorms(tp.time,Vectors(1,phasefield),gNorm))
      return false;
    else if (!gNorm.empty())
    {
      norm = gNorm.front();
      norm.push_back(norm(2));
      if (tp.step == 1)
        eps_d0 = norm(2);
      norm(2) -= eps_d0;
      if (norm.size() > 0 && utl::trunc(norm(1)) != 0.0)
        IFEM::cout <<"  L2-norm: |c^h| = (c^h,c^h)^0.5 : "
                   << sqrt(norm(1)) << std::endl;
      if (norm.size() > 2 && utl::trunc(norm(3)) != 0.0)
        IFEM::cout <<"  Dissipated energy:       eps_d : "
                   << norm(3) << std::endl;
    }

    return true;
  }

  //! \brief Sets initial conditions.
  void setInitialConditions() { SIM::setInitialConditions(*this); }

  //! \brief Sets the tensile energy vector from the elasticity problem.
  void setTensileEnergy(const RealArray* te)
  {
    static_cast<CahnHilliard*>(Dim::myProblem)->setTensileEnergy(te);
  }

  //! \brief Returns a const reference to the global norms.
  const Vector& getGlobalNorms() const { return norm; }

  //! \brief Returns a const reference to the current solution.
  const Vector& getSolution(int = 0) { return phasefield; }

  //! \brief Returns the maximum number of iterations (unlimited).
  int getMaxit() const { return 9999; }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //!
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    return this->solveStep(tp,false) ? SIM::CONVERGED : SIM::DIVERGED;
  }

protected:
  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"cahnhilliard"))
      return this->Dim::parse(elem);

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
        Dim::myProblem->parse(child);
        this->Dim::parse(child);
      }

    return true;
  }

private:
  Vector phasefield; //!< Current phase field solution
  Matrix projSol;    //!< Projected solution fields
  Vector norm;       //!< Global norm values
  double eps_d0;     //!< Initial eps_d value, subtracted from following values
};

#endif
