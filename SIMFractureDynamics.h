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
#include <fstream>


/*!
  \brief Driver class for fracture dynamics simulators.
  \details A fracture dynamics simulator is a coupling between
  a dynamic elasticity solver and a phase field solver.
*/

template<class SolidSolver, class PhaseSolver>
class SIMFracture : public SIMCoupled<SolidSolver,PhaseSolver>
{
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMFracture(SolidSolver& s1, PhaseSolver& s2)
    : SIMCoupled<SolidSolver,PhaseSolver>(s1,s2) {}
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

  //! \brief Saves the converged results to VTF-file of a given time step.
  //! \details It also writes global energy quantities to file for plotting.
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (!energFile.empty() && this->S1.getProcessAdm().getProcId() == 0)
    {
      std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);

      if (tp.step == 1)
        os <<"#t eps_e external_energy eps+ eps- eps_b |c| eps_d"<< std::endl;

      const Vector& n1 = this->S1.getGlobalNorms();
      const Vector& n2 = this->S2.getGlobalNorms();

      size_t i;
      os << std::setprecision(11) << std::setw(6) << std::scientific
         << tp.time.t;
      for (i = 0; i < n1.size(); i++) os <<" "<< n1[i];
      for (i = 0; i < n2.size(); i++) os <<" "<< n2[i];
      os << std::endl;
    }

    return this->S2.saveStep(tp,nBlock) && this->S1.saveStep(tp,nBlock);
  }

  //! \brief Assigns the file name for global energy output.
  void setEnergyFile(const char* fName)
  {
    if (fName)
    {
      energFile = fName;
      IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
    }
  }

private:
  std::string energFile; //!< File name for global energy output
};

#endif
