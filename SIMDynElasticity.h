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

#include "NonLinSIM.h"


class Elasticity;
class RealFunc;


/*!
  \brief Linear quasi-static solution driver.
*/

class LinSIM : public NonLinSIM
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit LinSIM(SIMbase& sim) : NonLinSIM(sim,NonLinSIM::NONE) {}
  //! \brief Empty destructor.
  virtual ~LinSIM() {}
};



/*!
  \brief Driver class for dynamic elasticity problems with fracture.
*/

template< class Dim, class DynSIM, class Sim>
class SIMDynElasticity : public Sim
{
public:
  //! \brief Default constructor.
  SIMDynElasticity();

  //! \brief Constructor for mixed problems.
  explicit SIMDynElasticity(const std::vector<unsigned char>& nf);

  //! \brief Empty destructor.
  virtual ~SIMDynElasticity() {}

  //! \brief Prints out problem-specific data to the log stream.
  void printProblem() const override;

  //! \brief Initializes the problem.
  bool init(const TimeStep& tp, bool = false) override;

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running result block counter
  bool saveStep(const TimeStep& tp, int& nBlock) override;

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep& tp) override { return dSim.advanceStep(tp,false); }

  //! \brief Prints out time step identification.
  //! \param[in] istep Time step counter
  //! \param[in] time Parameters for time-dependent simulations
  void printStep(int istep, const TimeDomain& time) const override;

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp) override;

  //! \brief Computes solution norms, etc. on the converged solution.
  bool postSolve(TimeStep& tp);

  //! \brief Updates the strain energy density for the current solution.
  bool updateStrainEnergyDensity(const TimeStep& tp);

  //! \brief Returns the tensile energy in gauss points.
  const RealArray* getTensileEnergy() const;

  //! \brief Returns a const reference to the global norms.
  const Vector& getGlobalNorms() const { return gNorm; }

  //! \brief Dummy method.
  void parseStaggering(const TiXmlElement*) {}

  //! \brief Assigns the file name for global energy output.
  void setEnergyFile(const char* fName);

  //! \brief Returns a const reference to current solution vector.
  const Vector& getSolution(int i) const override { return dSim.getSolution(i);}
  //! \brief Returns a const reference to the solution vectors.
  const Vectors& getSolutions() const override { return dSim.getSolutions(); }
  //! \brief Returns a reference to the solution vectors (for assignment).
  Vectors& theSolutions() override { return dSim.theSolutions(); }

  //! \brief Updates the solution vectors.
  void setSolutions(const Vectors& dvec);

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //! \param[in] stage Option, 1: solve only the first iteration,
  //! 2: solve for the remaining iterations, else: solve the whole time step
  SIM::ConvStatus solveIteration(TimeStep& tp, char stage = 0);

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return dSim.getMaxit(); }

  //! \brief Checks whether an internal crack pressure has been specified.
  RealFunc* haveCrackPressure() const;

  //! \brief Serializes current internal state for restarting purposes.
  bool serialize(SIMsolution::SerializeMap& data) const override;

  //! \brief Restores the internal state from serialized data.
  bool deSerialize(const SIMsolution::SerializeMap& data) override;

protected:
  //! \brief Returns the actual integrand.
  Elasticity* getIntegrand() override;

  using Sim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override;

private:
  std::string energFile; //!< File name for global energy output

  DynSIM dSim;    //!< Dynamic solution driver
  Matrix projSol; //!< Projected secondary solution fields
  Matrix eNorm;   //!< Element norm values
  Vector gNorm;   //!< Global norm values
  int    vtfStep; //!< VTF file step counter
  int    outPrec; //!< Precision on solution norm outputs
};

#endif
