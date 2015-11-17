// $Id$
//==============================================================================
//!
//! \file SIMDynElasticity.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Wrapper equipping the elasticity solver with
//! time-stepping support and phase field coupling.
//!
//==============================================================================

#ifndef _SIM_DYN_ELASTICITY_H_
#define _SIM_DYN_ELASTICITY_H_

#include "NewmarkSIM.h"
#include "SIMElasticity.h"
#include "FractureElasticity.h"
#include "DataExporter.h"


/*!
  \brief Driver wrapping an elasticity solver with an ISolver interface.
*/

template<class Dim, class DynSIM=NewmarkSIM>
class SIMDynElasticity : public SIMElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMDynElasticity() : SIMElasticity<Dim>(false),
                       fracEl(Dim::dimension), dSim(*this)
  {
    Dim::myHeading = "Elasticity solver";
    Dim::myProblem = &fracEl;
  }

  //! \brief The destructor nullifies the integrand pointers.
  virtual ~SIMDynElasticity()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

  //! \brief Prints out problem-specific data to the log stream.
  virtual void printProblem() const
  {
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      dSim.printProblem();
    else
      this->SIMElasticity<Dim>::printProblem();
    --ncall;
  }

  //! \brief Registers fields for data output.
  void registerFields(DataExporter& exporter)
  {
    exporter.registerField("u","solid displacement",DataExporter::SIM,
                           DataExporter::PRIMARY | DataExporter::SECONDARY);
    exporter.setFieldValue("u",this,&dSim.getSolution());
  }

  //! \brief Initializes the problem.
  bool init(const TimeStep&)
  {
    dSim.initPrm();
    dSim.initSol(3);

    this->setMode(SIM::INIT);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
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
    return dSim.saveModel(geoBlk,nBlock,fileName);
  }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    double old = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    bool ok = this->savePoints(dSim.getSolution(),tp.time.t,tp.step);
    utl::zero_print_tol = old;

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0 || !ok)
      return ok;

    return dSim.saveStep(tp.step/Dim::opt.saveInc,nBlock,tp.time.t);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp) { return dSim.advanceStep(tp,false); }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (Dim::msgLevel >= 1)
      IFEM::cout <<"\n  Solving the elasto-dynamics problem...";
    return dSim.solveStep(tp) == SIM::CONVERGED;
  }

  //! \brief Returns the tensile energy in gauss points.
  const double* getTensileEnergy() const { return fracEl.getTensileEnergy(); }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand() { return &fracEl; }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    bool result = true;
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      result = dSim.parse(elem);
    else
      result = this->SIMElasticity<Dim>::parse(elem);
    --ncall;
    return result;
  }

private:
  FractureElasticity fracEl; //!< Problem definition
  DynSIM             dSim;   //!< Dynamic solution driver
};

#endif
