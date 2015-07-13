// $Id$
//==============================================================================
//!
//! \file SIMElasticityWrap.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Wrapper equpping the nonlinear elasticity solver with
//! time-stepping support and phase field coupling.
//!
//==============================================================================

#ifndef _SIM_ELASTICITY_WRAP_H_
#define _SIM_ELASTICITY_WRAP_H_

#include "SIMSolver.h"
#include "DataExporter.h"
#include "Profiler.h"
#include "TimeStep.h"
#include "ASMstruct.h"
#include "HHTSIM.h"
#include "SIMFiniteDefEl.h"


/*!
  \brief Driver wrapping a nonlinear elasticity solver with an ISolver interface.
*/

template<class Dim, class NLSIM=HHTSIM> class SIMElasticityWrap : public SIMFiniteDefEl<Dim>
{
public:
  typedef bool SetupProps;

  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  SIMElasticityWrap() : SIMFiniteDefEl<Dim>(false, std::vector<int>()), nSim(*this)
  {
    this->myHeading = "Elasticity solver";
  }
  //! \brief Destructor.
  virtual ~SIMElasticityWrap() { this->setVTF(NULL); }

  //! \brief Registers fields for output to a data exporter.
  virtual void registerFields(DataExporter& exporter)
  {
    exporter.registerField("solid displacement", "solid displacement",
                           DataExporter::SIM, DataExporter::PRIMARY |
                                              DataExporter::SECONDARY);
    exporter.setFieldValue("solid displacement", this, &sol);
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  virtual bool saveModel(char* fileName, int& geoBlk, int& nBlock) { return true; }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    double old = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    this->savePoints(sol,tp.time.t,tp.step);
    utl::zero_print_tol = old;

    if (Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    return this->writeGlvS(sol,iDump,nBlock,tp.time.t);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp) { return nSim.advanceStep(tp, true); }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    return nSim.solveStep(tp);
  }

  bool postSolve(const TimeStep&, bool) { return true; }

  bool init(const TimeStep&)
  {
    sol.resize(this->getNoDOFs(), true);
    return true;
  }

  const Vector* getTensileEnergy() const { return nullptr; }
protected:
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"elasticity"))
      return nSim.parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      this->getIntegrand()->parse(child);

    return nSim.parse(elem);
  }

  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    return nullptr;
  }

private:
  Vector sol;
  NLSIM nSim;
};


//! \brief Partial specialization for configurator
  template<class Dim>
struct SolverConfigurator< SIMElasticityWrap<Dim> > {
  int setup(SIMElasticityWrap<Dim>& ad,
            const typename SIMElasticityWrap<Dim>::SetupProps& props, char* infile)
  {
    utl::profiler->start("Model input");

    // Reset the global element and node numbers
    ASMstruct::resetNumbering();
    if (!ad.read(infile))
      return 2;

    utl::profiler->stop("Model input");

    // Preprocess the model and establish data structures for the algebraic system
    if (!ad.preprocess())
      return 3;

    // Initialize the linear solvers
    ad.setMode(SIM::STATIC);
    ad.initSystem(ad.opt.solver,1,1);
    ad.setAssociatedRHS(0,0);
    ad.setQuadratureRule(ad.opt.nGauss[0]);
    ad.init(TimeStep());

    return 0;
  }
};

#endif
