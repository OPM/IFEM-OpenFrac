// $Id$
//==============================================================================
//!
//! \file SIMExplPhaseField.h
//!
//! \date May 27 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver representing an explicit phase-field.
//!
//==============================================================================

#ifndef _SIM_EXPL_PHASE_FIELD_H
#define _SIM_EXPL_PHASE_FIELD_H

#include "SIMbase.h"
#include "SIMdummy.h"
#include "SIMenums.h"

class SIMoutput;
class DataExporter;
class TimeStep;
class VTF;
namespace LR { class LRSpline; class RefineData; }


/*!
  \brief Driver class for an explicit phase-field.
*/

class SIMExplPhaseField : public SIMdummy<SIMbase>
{
public:
  //! \brief Default constructor.
  explicit SIMExplPhaseField(SIMoutput* gridOwner = nullptr);
  //! \brief The destructor deletes the explicit phase field function.
  virtual ~SIMExplPhaseField();

  //! \brief Registers fields for data output.
  void registerFields(DataExporter& exporter);

  //! \brief Initializes the problem.
  bool init(const TimeStep&);

  //! \brief Saves the converged results of a given time step to VTF file.
  bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool = true);

  //! \brief Returns the initial crack function.
  RealFunc* getInitCrack() const { return phaseFunc; }

  //! \brief Dummy method.
  bool postSolve(TimeStep&) { return true; }
  //! \brief Dummy method.
  bool advanceStep(TimeStep&) { return true; }
  //! \brief Dummy method.
  bool serialize(std::map<std::string,std::string>&) { return false; }
  //! \brief Dummy method.
  bool deSerialize(const std::map<std::string,std::string>&) { return false; }
  //! \brief Dummy method.
  bool dumpGeometry(std::ostream& os) const { return false; }
  //! \brief Dummy method.
  bool saveResidual(const TimeStep&, const Vector&, int&) { return true; }
  //! \brief Dummy method.
  bool checkStopCriterion () const { return false; }
  //! \brief Dummy method.
  void setTensileEnergy(const RealArray*) {}
  //! \brief Dummy method.
  void setVTF(VTF*) {}
  //! \brief Dummy method.
  double getNorm(Vector&, size_t = 0) const { return 0.0; }
  //! \brief Dummy method.
  const Vector& getGlobalNorms() const { static Vector g; return g; }
  //! \brief Returns a const reference to the current solution.
  const Vector& getSolution(int = 0) const { return phaseField; }
  //! \brief Updates the solution vector.
  void setSolution(const Vector& vec) { phaseField = vec; }
  //! \brief Returns the maximum number of iterations (unlimited).
  int getMaxit() const { return 9999; }
  //! \brief Dummy method.
  int getOutPrec() const { return 0; }
  //! \brief Dummy method.
  SIM::ConvStatus solveIteration(TimeStep&) { return SIM::CONVERGED; }
  //! \brief Dummy method.
  Vector getHistoryField() const { return Vector(); }
  //! \brief Dummy method.
  void setHistoryField(const RealArray&) {}
  //! \brief Dummy method.
  bool refine(const LR::RefineData&) { return false; }
  //! \brief Dummy method.
  void getBasis(std::vector<LR::LRSpline*>&) {}
  //! \brief Dummy method.
  bool transferHistory(const RealArray&,
                       std::vector<LR::LRSpline*>&) { return true; }

protected:
  using SIMbase::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem);

private:
  SIMoutput* myOwner;    //!< The FE mesh holder
  int        myStep;     //!< VTF file step counter
  RealFunc*  phaseFunc;  //!< Explicit phase field function
  Vector     phaseField; //!< Nodal phase field values
};

#endif
