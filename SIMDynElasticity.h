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
#include "FractureElasticityVoigt.h"
#include "DataExporter.h"


/*!
  \brief Driver wrapping an elasticity solver with an ISolver interface.
*/

template<class Dim, class DynSIM=NewmarkSIM>
class SIMDynElasticity : public SIMElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMDynElasticity(bool voigt = false) : SIMElasticity<Dim>(false), dSim(*this)
  {
    Dim::myHeading = "Elasticity solver";
    if (voigt)
      Dim::myProblem = new FractureElasticityVoigt(Dim::dimension);
    else
      Dim::myProblem = new FractureElasticity(Dim::dimension);
  }

  //! \brief Empty destructor.
  virtual ~SIMDynElasticity() {}

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

    int iDump = tp.step/Dim::opt.saveInc;
    if (!dSim.saveStep(iDump,nBlock,tp.time.t))
      return false;

    // Write projected solution fields to VTF-file
    if (!Dim::opt.project.empty())
      if (!this->writeGlvP(projSol,iDump,nBlock,110,
                           Dim::opt.project.begin()->second.c_str()))
        return false;

    // Write element norms to VTF-file
    return this->writeGlvN(eNorm,iDump,nBlock);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp) { return dSim.advanceStep(tp,false); }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (Dim::msgLevel >= 1)
      IFEM::cout <<"\n  Solving the elasto-dynamics problem...";

    if (dSim.solveStep(tp) != SIM::CONVERGED)
      return false;

    // Update strain energy density for the converged solution
    this->setMode(SIM::RECOVERY);
    if (!this->assembleSystem(tp.time,dSim.getSolutions()))
      return false;

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
      const Vector& gNorm = gNorms.front();
      if (gNorm.size() > 0 && utl::trunc(gNorm(1)) != 0.0)
        IFEM::cout <<"  Elastic strain energy:           eps_e : "
                   << gNorm(1) << std::endl;
      if (gNorm.size() > 1 && utl::trunc(gNorm(2)) != 0.0)
        IFEM::cout <<"  External energy: ((f,u^h)+(t,u^h))^0.5 : "
                   << sqrt(gNorm(2)) << std::endl;
    }

    return true;
  }

  //! \brief Returns the tensile energy in gauss points.
  const RealArray* getTensileEnergy() const
  {
    return static_cast<FractureElasticity*>(Dim::myProblem)->getTensileEnergy();
  }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    return static_cast<Elasticity*>(Dim::myProblem);
  }

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
  DynSIM dSim;    //!< Dynamic solution driver
  Matrix projSol; //!< Projected secondary solution fields
  Matrix eNorm;   //!< Element norm values
};

#endif
