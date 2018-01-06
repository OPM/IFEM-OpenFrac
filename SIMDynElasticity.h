// $Id$
//==============================================================================
//!
//! \file SIMDynElasticity.h
//!
//! \date Dec 04 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dynamic simulation driver for elasticity problems with fracture.
//!
//==============================================================================

#ifndef _SIM_DYN_ELASTICITY_H_
#define _SIM_DYN_ELASTICITY_H_

#include "NewmarkSIM.h"
#include "SIMElasticity.h"
#include "FractureElasticityVoigt.h"
#include "DataExporter.h"
#include <fstream>


/*!
  \brief Driver class for dynamic elasticity problems with fracture.
*/

template<class Dim, class DynSIM=NewmarkSIM>
class SIMDynElasticity : public SIMElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMDynElasticity() : dSim(*this)
  {
    Dim::myHeading = "Elasticity solver";
    vtfStep = 0;
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
    int flag = DataExporter::PRIMARY;
    if (!Dim::opt.pSolOnly)
      flag |= DataExporter::SECONDARY;
    exporter.registerField("u","solid displacement",DataExporter::SIM,flag);
    exporter.setFieldValue("u",this,&dSim.getSolution());
  }

  //! \brief Initializes the problem.
  bool init(const TimeStep&)
  {
    dSim.initPrm();
    dSim.initSol(dynamic_cast<NewmarkSIM*>(&dSim) ? 3 : 1);

    bool ok = this->setMode(SIM::INIT);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    this->registerField("solution",dSim.getSolution());
    return this->setInitialConditions() && ok;
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    return dSim.saveModel(geoBlk,nBlock,fileName);
  }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running result block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    double old = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    bool ok = this->savePoints(dSim.getSolution(),tp.time.t,tp.step);
    utl::zero_print_tol = old;

    if (!energFile.empty() && tp.step > 0 && Dim::adm.getProcId() == 0)
    {
      std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);

      Vector Bforce, Rforce;
      this->getBoundaryForce(Bforce,dSim.getSolutions(),tp);
      this->getBoundaryReactions(Rforce);

      if (tp.step == 1)
      {
        size_t i;
        os <<"#t eps_e external_energy eps+ eps- eps_b";
        for (i = 0; i < Bforce.size(); i++)
          os <<" load_"<< char('X'+i);
        for (i = 0; i < Rforce.size(); i++)
          os <<" react_"<< char('X'+i);
        os << std::endl;
      }

      os << std::setprecision(11) << std::setw(6) << std::scientific
         << tp.time.t;
      for (double n : gNorm)  os <<" "<< n;
      for (double f : Bforce) os <<" "<< utl::trunc(f);
      for (double f : Rforce) os <<" "<< utl::trunc(f);
      os << std::endl;
    }

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0 || !ok)
      return ok;

    // Write primary and secondary (of requested) solution fields to VTF-file
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
  virtual bool advanceStep(TimeStep& tp) { return dSim.advanceStep(tp,false); }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
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
    RealArray RF;
    if (this->getCurrentReactions(RF,dSim.getSolution()))
    {
      IFEM::cout <<"  Total reaction forces:          Sum(R) :";
      for (size_t i = 1; i < RF.size(); i++)
        IFEM::cout <<" "<< utl::trunc(RF[i]);
      double Ru = RF.front();
      if (utl::trunc(Ru) != 0.0)
        IFEM::cout <<"\n  displacement*reactions:          (R,u) : "<< Ru;
      IFEM::cout << std::endl;
    }

    // Update strain energy density for the converged solution
    if (!this->updateStrainEnergyDensity(tp))
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
                   << (gNorm(2) < 0.0 ? -sqrt(-gNorm(2)) : sqrt(gNorm(2)))
                   << std::endl;
    }

    return true;
  }

  //! \brief Updates the strain energy density for the current solution.
  bool updateStrainEnergyDensity(const TimeStep& tp)
  {
    this->setMode(SIM::RECOVERY);
    return this->assembleSystem(tp.time,dSim.getSolutions());
  }

  //! \brief Returns the tensile energy in gauss points.
  const RealArray* getTensileEnergy() const
  {
    return static_cast<Elasticity*>(Dim::myProblem)->getTensileEnergy();
  }

  //! \brief Returns a const reference to the global norms.
  const Vector& getGlobalNorms() const { return gNorm; }

  //! \brief Dummy method.
  void parseStaggering(const TiXmlElement*) {}

  //! \brief Assigns the file name for global energy output.
  void setEnergyFile(const char* fName)
  {
    if (fName)
    {
      energFile = fName;
      IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
    }
  }

  //! \brief Returns a const reference to current solution vector.
  const Vector& getSolution(int idx = 0) const { return dSim.getSolution(idx); }
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
  //! \param[in] stage Option, 1: solve only the first iteration,
  //! 2: solve for the remaining iterations, else: solve the whole time step
  SIM::ConvStatus solveIteration(TimeStep& tp, char stage = 0)
  {
    if (Dim::msgLevel == 1 && tp.iter == 0 && stage == 1)
      IFEM::cout <<"\n  step="<< tp.step <<"  time="<< tp.time.t << std::endl;
    dSim.setSubIteration(tp.iter == 0 ? DynSIM::FIRST : DynSIM::ITER);

    if (tp.iter == 0 && stage == 1)
      return dSim.solveIteration(tp);
    else if (tp.iter > 0 && stage == 2)
    {
      SIM::ConvStatus status = SIM::OK;
      while (tp.iter <= dSim.getMaxit())
        switch ((status = dSim.solveIteration(tp))) {
        case SIM::OK:
        case SIM::SLOW:
          tp.iter++;
          break;
        default:
          return status;
	}
      return SIM::DIVERGED; // No convergence in maxit iterations
    }

    // Solve the whole time step
    TimeStep myTp(tp); // Make a copy to avoid destroying the iteration counter
    return dSim.solveStep(myTp);
  }

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return dSim.getMaxit(); }

  //! \brief Checks whether an internal crack pressure has been specified.
  bool haveCrackPressure() const
  {
    FractureElasticity* fel = dynamic_cast<FractureElasticity*>(Dim::myProblem);
    return fel ? fel->getCrackPressure() != nullptr : false;
  }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem) // Using the Voigt formulation by default
      Dim::myProblem = new FractureElasticityVoigt(Dim::dimension);
    return static_cast<Elasticity*>(Dim::myProblem);
  }

  using SIMElasticity<Dim>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    bool result = true;
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      result = dSim.parse(elem);
    else if (!strcasecmp(elem->Value(),"elasticity") && !Dim::myProblem)
    {
      std::string form("voigt");
      if (utl::getAttribute(elem,"formulation",form,true) && form != "voigt")
        Dim::myProblem = new FractureElasticity(Dim::dimension);
      result = this->SIMElasticity<Dim>::parse(elem);

      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        if (!strcasecmp(child->Value(),"noGeometricStiffness"))
          this->setIntegrationPrm(3,1); // Disable geometric stiffness
    }
    else
      result = this->SIMElasticity<Dim>::parse(elem);
    --ncall;
    return result;
  }

private:
  std::string energFile; //!< File name for global energy output

  DynSIM dSim;    //!< Dynamic solution driver
  Matrix projSol; //!< Projected secondary solution fields
  Matrix eNorm;   //!< Element norm values
  Vector gNorm;   //!< Global norm values
  int    vtfStep; //!< VTF file step counter
};

#endif
