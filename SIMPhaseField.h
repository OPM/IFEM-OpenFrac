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
  SIMPhaseField(char order = 2) : Dim(1)
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
  bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMPhaseField::solveStep");

    if (Dim::msgLevel == 1)
      IFEM::cout <<"\n  Solving crack phase field at step="<< tp.step
                 <<" time="<< tp.time.t << std::endl;

    this->setMode(SIM::STATIC);
    if (!this->assembleSystem())
      return false;

    if (!this->solveSystem(phasefield,1,
                           Dim::msgLevel > 1 ? "phasefield  " : nullptr))
      return false;

    static_cast<CahnHilliard*>(Dim::myProblem)->clearInitialCrack();
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
      if (norm.size() > 0 && utl::trunc(norm(1)) != 0.0)
        IFEM::cout <<"  L2-norm: |c^h| = (c^h,c^h)^0.5 : "
                   << sqrt(norm(1)) << std::endl;
      if (norm.size() > 1 && utl::trunc(norm(2)) != 0.0)
        IFEM::cout <<"  Dissipated energy:       eps_d : "
                   << norm(2) << std::endl;
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
};

#endif
