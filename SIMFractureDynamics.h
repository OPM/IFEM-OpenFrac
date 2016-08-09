// $Id$
//==============================================================================
//!
//! \file SIMFactureDynamics.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for fracture-dynamic problems.
//!
//==============================================================================

#ifndef _SIM_FRACTURE_DYNAMICS_H_
#define _SIM_FRACTURE_DYNAMICS_H_

#include "SIMCoupled.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "LRSpline/LRSplineSurface.h"
#endif
#include <fstream>
#include <numeric>


/*!
  \brief Driver class for fracture dynamics simulators.
  \details A fracture dynamics simulator is a coupling between
  a dynamic elasticity solver and a phase field solver.
*/

template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling=SIMCoupled>
class SIMFracture : public Coupling<SolidSolver,PhaseSolver>
{
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMFracture(SolidSolver& s1, PhaseSolver& s2, const std::string& inputfile)
    : Coupling<SolidSolver,PhaseSolver>(s1,s2), infile(inputfile), aMin(0.0) {}
  //! \brief Empty destructor.
  virtual ~SIMFracture() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies()
  {
    this->S1.registerDependency(&this->S2,"phasefield",1);
    // The tensile energy is defined on integration points and not nodal points.
    // It is a global buffer array across all patches in the model.
    // Use an explicit call instead of normal couplings for this.
    this->S2.setTensileEnergy(this->S1.getTensileEnergy());
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    if (tp.step == 1 && this->S1.haveCrackPressure())
      // Start the initial step by solving the phase-field first
      if (!this->S2.solveStep(tp,false))
        return false;

    return this->Coupling<SolidSolver,PhaseSolver>::solveStep(tp);
  }

  //! \brief Saves the converged results to VTF-file of a given time step.
  //! \details It also writes global energy quantities to file for plotting.
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (!energFile.empty() && this->S1.getProcessAdm().getProcId() == 0)
    {
      std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);

      size_t i;
      Vector BF, RF;
      this->S1.getBoundaryForce(BF,this->S1.getSolutions(),tp);
      this->S1.getBoundaryReactions(RF);

      if (tp.step == 1)
      {
        os <<"#t eps_e external_energy eps+ eps- eps_b |c|"
           <<" eps_d-eps_d(0) eps_d";
        for (i = 0; i < BF.size(); i++)
          os <<" load_"<< char('X'+i);
        for (i = 0; i < RF.size(); i++)
          os <<" react_"<< char('X'+i);
        os << std::endl;
      }

      const Vector& n1 = this->S1.getGlobalNorms();
      const Vector& n2 = this->S2.getGlobalNorms();

      os << std::setprecision(11) << std::setw(6) << std::scientific
         << tp.time.t;
      for (size_t i = 0; i < n1.size(); i++) os <<" "<< n1[i];
      os <<" "<< (n2.size() > 2 ? n2[1] : 0.0);
      os <<" "<< (n2.size() > 1 ? n2[n2.size()-2] : 0.0);
      os <<" "<< (n2.size() > 0 ? n2.back() : 0.0);
      for (i = 0; i < BF.size(); i++)
        os <<" "<< utl::trunc(BF[i]);
      for (i = 0; i < RF.size(); i++)
        os <<" "<< utl::trunc(RF[i]);
      os << std::endl;
    }

    return this->S2.saveStep(tp,nBlock) && this->S1.saveStep(tp,nBlock);
  }

  //! \brief Parses sub-iteration parameters from an XML element.
  virtual void parseSubiteration(const TiXmlElement*) {}

  //! \brief Assigns the file name for global energy output.
  void setEnergyFile(const char* fName)
  {
    if (fName)
    {
      energFile = fName;
      IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
    }
  }

  //! \brief Stores current solution state in an internal buffer.
  void saveState()
  {
    sols = this->S1.getSolutions();
    sols.push_back(this->S2.getSolution());
    hsol = this->S2.getHistoryField();
  }

  //! \brief Refines the mesh on the initial configuration.
  bool initialRefine(double beta, double min_frac, int nrefinements)
  {
    if (this->S2.getInitRefine() >= nrefinements)
      return true; // Grid is sufficiently refined during input parsing

    TimeStep step0;
    int newElements = 1;
    for (step0.iter = 0; newElements > 0; step0.iter++)
      if (!this->S2.solveStep(step0))
        return false;
      else
        newElements = this->adaptMesh(beta,min_frac,nrefinements);

    return newElements == 0;
  }

  //! \brief Refines the mesh with transfer of solution onto the new mesh.
  int adaptMesh(double beta, double min_frac, int nrefinements)
  {
#ifdef HAS_LRSPLINE
    ASMu2D* pch = dynamic_cast<ASMu2D*>(this->S1.getPatch(1));
    if (!pch)
      return -1;

    if (aMin <= 0.0) // maximum refinements per element
    {
      double redMax = pow(2.0,nrefinements);
      aMin = pch->getBasis()->getElement(0)->area()/(redMax*redMax);
    }

    // Fetch element norms to use as refinement criteria
    Vector eNorm;
    double gNorm = this->S2.getNorm(eNorm,3);

    // Sort element indices based on comparing values in eNorm
    IntVec idx(eNorm.size());
    std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),
              [&eNorm](size_t i1, size_t i2) { return eNorm[i1] < eNorm[i2]; });

    double eMin = min_frac < 0.0 ? -min_frac*gNorm/sqrt(idx.size()) : min_frac;
    size_t eMax = beta < 0.0 ? idx.size() : idx.size()*beta/100.0;
    IFEM::cout <<"\n  Lowest element: "<< std::setw(8) << idx.front()
               <<"    |c| = "<< eNorm[idx.front()]
               <<"\n  Highest element:"<< std::setw(8) << idx.back()
               <<"    |c| = "<< eNorm[idx.back()]
               <<"\n  Minimum |c|-value for refinement: "<< eMin
               <<"\n  Minimum element area: "<< aMin << std::endl;

    IntVec elements; // Find the elements to refine
    for (size_t i = 0; i < idx.size() && elements.size() < eMax; i++)
      if (eNorm[idx[i]] > eMin)
        break;
      else if (pch->getBasis()->getElement(idx[i])->area() > aMin+1.0e-12)
        elements.push_back(idx[i]);

    if (elements.empty())
      return 0;

    IFEM::cout <<"  Elements to refine: "<< elements.size()
               <<" (|c| = ["<< eNorm[elements.front()]
               <<","<< eNorm[elements.back()] <<"])\n"<< std::endl;

    LR::LRSplineSurface* oldBasis = nullptr;
    if (!hsol.empty()) oldBasis = pch->getBasis()->copy();

    // Do the mesh refinement
    LR::RefineData prm;
    prm.options = { 10, 1, 2, 0, 1 };
    prm.elements = pch->getFunctionsForElements(elements);
    if (!this->S1.refine(prm,sols) || !this->S2.refine(prm))
      return -2;

    // Re-initialize the simulators for the new mesh
    this->S1.clearProperties();
    this->S2.clearProperties();
    if (!this->S1.read(infile.c_str()) || !this->S2.read(infile.c_str()))
      return -3;

    if (!this->preprocess())
      return -4;

    if (!this->init(TimeStep()))
      return -5;

    if (!this->S1.initSystem(this->S1.opt.solver) ||
        !this->S2.initSystem(this->S2.opt.solver,1,1,false))
      return -6;

    // Transfer solution variables onto the new mesh
    if (!sols.empty())
    {
      IFEM::cout <<"\nTransferring "<< sols.size()-1 <<"x"<< sols.front().size()
                 <<" solution variables to new mesh for "<< this->S1.getName();
      this->S1.setSolutions(sols);
      IFEM::cout <<"\nTransferring "<< sols.back().size()
                 <<" solution variables to new mesh for "<< this->S2.getName();
      this->S2.setSolution(sols.back());
    }
    if (!hsol.empty())
    {
      IFEM::cout <<"\nTransferring "<< hsol.size()
                 <<" history variables to new mesh for "<< this->S2.getName()
                 << std::endl;
      this->S2.transferHistory2D(hsol,oldBasis);
      delete oldBasis;
    }

    return elements.size();
#else
    std::cerr <<" *** SIMFractureDynamics:adaptMesh: No LR-spline support.\n";
    return -1;
#endif
  }

private:
  std::string energFile; //!< File name for global energy output
  std::string infile;    //!< Input file parsed

  double    aMin; //!< Minimum element area
  Vectors   sols; //!< Solution state to transfer onto refined mesh
  RealArray hsol; //!< History field to transfer onto refined mesh
};


/*!
  \brief Driver class for quasi-static fracture simulators with sub-iterations.
*/

template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling=SIMCoupled>
class SIMFractureQstatic : public SIMFracture<SolidSolver,PhaseSolver,Coupling>
{
public:
  //! \brief The constructor forwards to the parent class contructor.
  SIMFractureQstatic(SolidSolver& s1, PhaseSolver& s2, const std::string& input)
    : SIMFracture<SolidSolver,PhaseSolver,Coupling>(s1,s2,input)
  {
    subIType = 1;
    maxSubIt = 50;
    subItTol = 1.0e-4;
  }

  //! \brief Empty destructor.
  virtual ~SIMFractureQstatic() {}

  //! \brief Parses sub-iteration parameters from an XML element.
  virtual void parseSubiteration(const TiXmlElement* elem)
  {
    this->S1.parseSubiteration(elem);
    utl::getAttribute(elem,"type",subIType);
    utl::getAttribute(elem,"tol",subItTol);
    if (subIType == 2)
      maxSubIt = this->S1.getMaxit();
    else
      utl::getAttribute(elem,"max",maxSubIt);
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    if (tp.step == 1 && this->S1.haveCrackPressure() && rHistory.empty())
      // Start the initial step by solving the phase-field first
      if (!this->S2.solveStep(tp,false))
        return false;

    return this->Coupling<SolidSolver,PhaseSolver>::solveStep(tp);
  }

  //! \brief Checks if the coupled simulator has converged.
  virtual SIM::ConvStatus checkConvergence(const TimeStep& tp,
                                           SIM::ConvStatus status1,
                                           SIM::ConvStatus status2)
  {
    if (status1 == SIM::FAILURE || status2 == SIM::FAILURE)
      return SIM::FAILURE;
    else if (status1 == SIM::DIVERGED || status2 == SIM::DIVERGED)
      return SIM::DIVERGED;
    else if (status1 != SIM::CONVERGED || status2 != SIM::CONVERGED)
      return SIM::OK;

    // Compute residual
    this->S1.setMode(SIM::RHS_ONLY);
    if (!this->S1.assembleSystem(tp.time,this->S1.getSolutions(),false))
      return SIM::FAILURE;

    Vector residual;
    if (!this->S1.extractLoadVec(residual))
      return SIM::FAILURE;

    int subIt = rHistory.size();
    double dummy, rNorm;
    this->S1.iterationNorms(Vector(), residual, dummy, rNorm, dummy);
    rHistory.push_back(rNorm);

    double beta = 10.0, rConv = rNorm/rHistory.front();
    IFEM::cout <<"  subit="<< subIt <<"  conv="<< rConv;
    if (subIt > 0)
    {
      double r0 = rHistory.front();
      double r2 = rHistory[subIt-1];
      beta = atan2(subIt*(r2-rNorm),r0-rNorm) * 180.0/M_PI;
      IFEM::cout <<"  beta="<< beta;
    }
    IFEM::cout << std::endl;

    if (rConv < subItTol)
    {
      rHistory.clear();
      return SIM::CONVERGED;
    }
    else if (subIt <= maxSubIt)
      return SIM::OK;

    std::cerr <<"SIMFractureQstatic::checkConvergence: Did not converge in "
              << maxSubIt <<" sub-iterations, bailing.."<< std::endl;
    return SIM::DIVERGED;
  }

private:
  int       subIType; //!< Sub-iteration type flag
  int       maxSubIt; //!< Maximum number of sub-iterations
  double    subItTol; //!< Sub-iteration tolerance
  RealArray rHistory; //!< Residual norm history for quasi-static simulations
};

#endif
