// $Id$
//==============================================================================
//!
//! \file SIMPhaseField.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for Cahn-Hilliard phase-field problems.
//!
//==============================================================================

#ifndef _SIM_PHASE_FIELD_H
#define _SIM_PHASE_FIELD_H

#include "SIMenums.h"
#include "SIMsolution.h"
#include "Vec3.h"


class CahnHilliard;
class DataExporter;
namespace LR { class LRSpline; }
class RealFunc;
class TimeStep;
class TiXmlElement;


/*!
  \brief Driver class for a Cahn-Hilliard phase-field simulator.
*/

template<class Dim> class SIMPhaseField : public Dim, public SIMsolution
{
  //! \brief Struct defining a user-specified plane for stopping the simulation.
  struct StopPlane
  {
    Vec3             normal; //!< Normal vector of the stop plane
    double           d;      //!< Distance offset of plane origin
    double           eps;    //!< Zero tolerance for nodes being on the plane
    std::vector<int> nodes;  //!< Nodal (control) points on the plane
    double           cstop;  //!< Stop when phase field value lower that this
  };

public:
  //! \brief Default constructor.
  SIMPhaseField(Dim* gridOwner = nullptr, size_t n = 2);

  //! \brief The destructor deletes the stop plane.
  virtual ~SIMPhaseField();

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "CahnHilliard"; }

  //! \brief Preprocessing performed before the FEM model generation.
  void preprocessA() override;

  //! \brief Preprocessing performed after the FEM model generation.
  bool preprocessB() override;

  //! \brief Registers fields for data output.
  void registerFields(DataExporter& exporter);

  //! \brief Initializes the problem.
  bool init(const TimeStep& tp);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock);

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Saves the force residual of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] residual Residual force vector
  //! \param nBlock Running VTF block counter
  bool saveResidual(const TimeStep& tp, const Vector& residual, int& nBlock);

  //! \brief Serializes current internal state for restarting purposes.
  bool serialize(SerializeMap& data) const override;

  //! \brief Restores the internal state from serialized data.
  bool deSerialize(const SerializeMap& data) override;

  //! \brief Advances the time step one step forward.
  bool advanceStep(const TimeStep&);

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool standalone = true);

  //! \brief Computes solution norms, etc. on the converged solution.
  bool postSolve(TimeStep& tp);

  //! \brief Prints a summary of the calculated solution to console.
  //! \param[in] solvec The solution vector
  //! \param[in] printSol Print solution only if size is less than this value
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] prec Number of digits after the decimal point in norm print
  void printSolutionSummary(const Vector& solvec, int printSol,
                            const char* compName, std::streamsize prec) override;

  //! \brief Returns \e true if terminating due to user-defined criteria.
  bool checkStopCriterion() const;

  //! \brief Returns the initial crack function.
  RealFunc* getInitCrack() const;

  //! \brief Sets the tensile energy vector from the elasticity problem.
  void setTensileEnergy(const RealArray* te);

  //! \brief Returns a list of element norm values.
  double getNorm(Vector& values, size_t idx = 1) const;

  //! \brief Returns a const reference to the global norms.
  const Vector& getGlobalNorms() const { return norm; }

  //! \brief Returns the maximum number of iterations (unlimited).
  int getMaxit() const { return 9999; }

  //! \brief Returns the norm output precision.
  int getOutPrec() const { return outPrec; }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //!
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp);

  //! \brief Returns the current history field.
  //! \details If projection has been done, the resulting control point values
  //! are returned, otherwise the Gauss point values are returned.
  RealArray getHistoryField() const;

  //! \brief Resets the history field from the provided array.
  void setHistoryField(const RealArray& hfield);

  //! \brief Extracts the LR-spline basis for the phase field
  void getBasis(std::vector<LR::LRSpline*>& basis);

  //! \brief Transfers history variables at Gauss/control points to new mesh.
  //! \param[in] oldH History variables associated with Gauss- or control points
  //! \param[in] oldBasis The LR-spline basis \a oldH is referring to
  bool transferHistory(const RealArray& oldH,
                       std::vector<LR::LRSpline*>& oldBasis);

  // Due to the multiple inheritance, the compiler needs to be told which
  // version of this method to use (even though they have different signature)
  using SIMsolution::getSolution;

protected:
  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override;

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override;

private:
  CahnHilliard* chp;  //!< The Cahn-Hilliard integrand
  StopPlane*    spln; //!< Stop simulation when crack penetrates this plane

  Matrix  projSol;    //!< Projected solution fields
  Matrix  eNorm;      //!< Element norm values
  Vector  norm;       //!< Global norm values
  double  eps_d0;     //!< Initial eps_d value, subtracted from following values
  int     vtfStep;    //!< VTF file step counter
  int     outPrec;    //!< Precision on solution norm outputs
  char    transferOp; //!< Solution transfer option
};

#endif
