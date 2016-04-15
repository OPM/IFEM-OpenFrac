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
#include "SIMPoroElasticity.h"
#include "FractureElasticityVoigt.h"


/*!
  \brief Driver wrapping an elasticity solver with an ISolver interface.
*/

template<class Dim, class DynSIM=NewmarkSIM>
class SIMDynElasticity : public SIMPoroElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMDynElasticity() : dSim(*this), vtfStep(0) {}

  //! \brief Constructor for mixed poroelastic problems.
  SIMDynElasticity(const std::vector<unsigned char>& flds)
    : SIMPoroElasticity<Dim>(flds), dSim(*this), vtfStep(0) {}

  //! \brief Empty destructor.
  virtual ~SIMDynElasticity() {}

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const
  {
    if (SIMElasticity<Dim>::myContext == "poroelasticity")
      return "PoroElasticity";
    else
      return "Elasticity";
  }

  //! \brief Prints out problem-specific data to the log stream.
  virtual void printProblem() const
  {
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      dSim.printProblem();
    else
      this->Dim::printProblem();
    --ncall;
  }

  //! \brief Initializes the solution vectors.
  virtual bool init(const TimeStep&)
  {
    dSim.initPrm();
    dSim.initSol(3);

    bool ok = this->setMode(SIM::INIT);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    return ok;
  }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running result block counter
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    double old = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    bool ok = this->savePoints(dSim.getSolution(),tp.time.t,tp.step);
    utl::zero_print_tol = old;

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0 || !ok)
      return ok;

    if (!dSim.saveStep(++vtfStep,nBlock,tp.time.t))
      return false;
    else if (tp.step < 1)
      return true;

    // Write projected solution fields to VTF-file
    if (!Dim::opt.project.empty())
      if (!this->writeGlvP(projSol,vtfStep,nBlock,110,
                           Dim::opt.project.begin()->second.c_str()))
        return false;

    // Write element norms to VTF-file
    return this->writeGlvN(eNorm,vtfStep,nBlock);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    dSim.advanceStep(tp,false);
    return this->SIMElasticity<Dim>::advanceStep(tp);
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    if (Dim::msgLevel >= 1)
      IFEM::cout <<"\n  Solving the elasto-dynamics problem...";

    if (dSim.solveStep(tp) != SIM::CONVERGED)
      return false;

    return this->postSolve(tp);
  }

  //! \brief Computes solution norms, etc. on the converged solution.
  bool postSolve(TimeStep& tp)
  {
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
      gNorm = gNorms.front();
      if (gNorm.size() > 0 && utl::trunc(gNorm(1)) != 0.0)
        IFEM::cout <<"  Elastic strain energy:           eps_e : "
                   << gNorm(1) << std::endl;
      if (gNorm.size() > 4 && utl::trunc(gNorm(5)) != 0.0)
        IFEM::cout <<"  Bulk energy:                     eps_b : "
                   << gNorm(5)
                   <<"\n  Tensile & compressive energies         : "
                   << gNorm(3) <<" "<< gNorm(4) << std::endl;
      if (gNorm.size() > 1 && utl::trunc(gNorm(2)) != 0.0)
        IFEM::cout <<"  External energy: ((f,u^h)+(t,u^h))^0.5 : "
                   << sqrt(gNorm(2)) << std::endl;
    }

    return true;
  }

  //! \brief Returns the tensile energy in gauss points.
  const RealArray* getTensileEnergy() const
  {
    return static_cast<Elasticity*>(Dim::myProblem)->getTensileEnergy();
  }

  //! \brief Returns a const reference to the global norms.
  const Vector& getGlobalNorms() const { return gNorm; }

  //! \brief Dummy method.
  void setEnergyFile(const std::string&) {}

  //! \brief Returns a const reference to current solution vector.
  virtual const Vector& getSolution(int i) const { return dSim.getSolution(i); }
  //! \brief Returns a const reference to the solution vectors.
  const Vectors& getSolutions() const { return dSim.getSolutions(); }

  //! \brief Updates the solution vectors.
  void setSolutions(const Vectors& dvec)
  {
    size_t nSol = dSim.getSolutions().size();
    for (size_t i = 0; i < nSol && i < dvec.size(); i++)
      dSim.setSolution(dvec[i],i);
  }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    return dSim.solveIteration(tp);
  }

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return dSim.getMaxit(); }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem) // Using the Voigt elasticity formulation by default
      Dim::myProblem = new FractureElasticityVoigt(Dim::dimension);
    return static_cast<Elasticity*>(Dim::myProblem);
  }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    bool result = true;
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      result = dSim.parse(elem);
    else if (!strcasecmp(elem->Value(),"elasticity"))
    {
      std::string form("voigt");
      if (!Dim::myProblem)
        if (utl::getAttribute(elem,"formulation",form,true) && form != "voigt")
          Dim::myProblem = new FractureElasticity(Dim::dimension);
      SIMElasticity<Dim>::myContext = "elasticity";
      result = this->SIMElasticity<Dim>::parse(elem);
    }
    else if (!strcasecmp(elem->Value(),"poroelasticity"))
    {
      if (!Dim::myProblem)
        Dim::myProblem = new PoroElasticity(Dim::dimension);
      SIMElasticity<Dim>::myContext = "poroelasticity";
      result = this->SIMElasticity<Dim>::parse(elem);
    }
    else
      result = this->Dim::parse(elem);

    --ncall;
    return result;
  }

private:
  DynSIM dSim;    //!< Dynamic solution driver
  Matrix projSol; //!< Projected secondary solution fields
  Matrix eNorm;   //!< Element norm values
  Vector gNorm;   //!< Global norm values
  int    vtfStep; //!< VTF file step counter
};

#endif
