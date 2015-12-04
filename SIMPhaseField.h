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

template<class Dim, class Integrand> class SIMPhaseField : public Dim
{
public:
  //! \brief Default constructor.
  SIMPhaseField() : Dim(1), cahnhill(Dim::dimension)
  {
    Dim::myHeading = "Cahn-Hilliard solver";
    Dim::myProblem = &cahnhill;
  }

  //! \brief The destructor nullifies the integrand pointers.
  virtual ~SIMPhaseField()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

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
    this->setMode(SIM::STATIC);
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

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (this->writeGlvS1(phasefield,iDump,nBlock,tp.time.t,"phasefield",9) < 0)
      return false;

    return this->writeGlvStep(iDump,tp.time.t);
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

    if (!this->assembleSystem())
      return false;

    if (!this->solveSystem(phasefield,1,
                           Dim::msgLevel > 1 ? "phasefield  " : NULL))
      return false;

    cahnhill.clearInitialCrack();

    return Dim::opt.project.empty() ? true : this->postSolve(tp);
  }

  //! \brief Projects the phase field onto the geometry basis of this simulator.
  bool postSolve(const TimeStep&, bool = false)
  {
    Matrix tmp;
    if (!this->project(tmp,phasefield,Dim::opt.project.begin()->first))
      return false;

    phasefield = tmp.getRow(1);
    return true;
  }

  //! \brief Sets initial conditions.
  void setInitialConditions() { SIM::setInitialConditions(*this); }

  //! \brief Sets the tensile energy vector from the elasticity problem.
  void setTensileEnergy(const double* te) { cahnhill.setTensileEnergy(te); }

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
  Integrand cahnhill;   //!< Problem definition
  Vector    phasefield; //!< Current phase field solution
};

#endif
