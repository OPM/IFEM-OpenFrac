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

#include "SIMElasticity.h"
#include "HHTSIM.h"
#include "DataExporter.h"
#include "Profiler.h"


/*!
  \brief Driver wrapping a nonlinear elasticity solver with an ISolver interface.
*/

template<class Dim, class NLSIM=HHTSIM>
class SIMElasticityWrap : public SIMElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMElasticityWrap() : SIMElasticity<Dim>(false), nSim(*this)
  {
    this->myHeading = "Elasticity solver";
  }
  //! \brief Empty destructor.
  virtual ~SIMElasticityWrap() {}

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
  bool saveModel(char* fileName, int& geoBlk, int& nBlock) { return true; }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
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
  virtual bool advanceStep(TimeStep& tp) { return nSim.advanceStep(tp,true); }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp) { return nSim.solveStep(tp); }

  //! \brief Initializes the solution vector.
  void initSol() { sol.resize(this->getNoDOFs(),true); }
  //! \brief Dummy method.
  bool init(const TimeStep&) { return true; }

  //! \brief Return the tensile energy in gauss points.
  const Vector* getTensileEnergy() const
  {
    return nullptr; // TODO
  } 

protected:
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),"elasticity"))
    {
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        this->getIntegrand()->parse(child);
    }
    return nSim.parse(elem);
  }

  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    return nullptr; // TODO
  }

private:
  Vector sol;
  NLSIM  nSim;
};

#endif
