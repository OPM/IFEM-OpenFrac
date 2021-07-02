// $Id$
//==============================================================================
//!
//! \file SIMFractureDynamics.h
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

#include "MatVec.h"

#include <string>


class RealFunc;
class SIMadmin;
class TimeStep;
class TiXmlElement;


/*!
  \brief Generic input file parser with profiling.
*/

bool readModel (SIMadmin& model, const std::string& infile);


/*!
  \brief Driver class for fracture dynamics simulators.
  \details A fracture dynamics simulator is a coupling between
  a dynamic elasticity solver and a phase field solver.
*/

template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
class SIMFracture : public Coupling<SolidSolver,PhaseSolver>
{
  //! Convenience type
  typedef Coupling<SolidSolver,PhaseSolver> CoupledSIM;

public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMFracture(SolidSolver& s1, PhaseSolver& s2, const std::string& inputfile);

  //! \brief Empty destructor.
  virtual ~SIMFracture() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies();

  //! \brief Returns \e true if the user-defined stop criterion has been met.
  bool stopped() const { return doStop; }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp);

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool firstS1 = true);

  //! \brief Saves the converged results to VTF-file of a given time step.
  //! \details It also writes global energy quantities to file for plotting.
  virtual bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Parses pre-refinement parameters from an XML element.
  void parsePreref(const TiXmlElement* elem);

  //! \brief Parses staggering parameters from an XML element.
  virtual void parseStaggering(const TiXmlElement* elem);

  //! \brief Assigns the file name for global energy output.
  void setEnergyFile(const char* fName);

  //! \brief Stores current solution state in an internal buffer.
  void saveState();

  //! \brief Restores the solution state from the internal buffer.
  void restoreState();

  //! \brief Refines the mesh based on an explicit mesh density function.
  bool preRefine(int nrefinements, int irefine, double refTol);

  //! \brief Refines the mesh on the initial configuration.
  bool initialRefine(double beta, double min_frac, int nrefinements);

  //! \brief Refines the mesh with transfer of solution onto the new mesh.
  int adaptMesh(double beta, double min_frac, int nrefinements,
                bool remeshOnly = false);

  //! \brief Dumps current mesh to the specified file.
  bool dumpMesh(const char* fileName);

protected:
  //! \brief Calculates and prints the solution and residual norms.
  double calcResidual(const TimeStep& tp, bool cycles = false);

public:
  //! \brief Returns solution vector to use for relaxation.
  virtual const Vector& getRelaxationVector() const
  { return this->S2.getSolution(); }
  //! \brief Updates the relaxed solution vector.
  virtual void setRelaxedSolution(const Vector& sol)
  { this->S2.setSolution(sol); }

  //! \brief Returns the residual of the elasticity equation.
  virtual const Vector& getAitkenResidual() const { return elastRes; }

private:
  std::string energFile; //!< File name for global energy output
  std::string infile;    //!< Input file parsed

  RealFunc* refFunc; //!< Pre-refinement function

  Vectors   sols; //!< Solution state to transfer onto refined mesh
  RealArray hsol; //!< History field to transfer onto refined mesh

  size_t irfStop; //!< Reaction force component to use as stop criterion
  double stopVal; //!< Stop simulation when less that this value

  double E0; //!< Energy norm of initial staggering cycle
  double Ec; //!< Energy norm of current staggering cycle
  double Ep; //!< Energy norm of previous staggering cycle

  Vector elastRes; //!< Residual force vector of the elasticity equation
  Vector residual; //!< Residual force vector of the phase field equation

protected:
  bool doStop; //!< If \e true, terminate due to user-defined stop criterion
};

#endif
