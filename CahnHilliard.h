// $Id$
//==============================================================================
//!
//! \file CahnHilliard.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Cahn Hilliard problems.
//!
//==============================================================================

#ifndef _CAHN_HILLIARD_H
#define _CAHN_HILLIARD_H

#include "IntegrandBase.h"


/*!
  \brief Class representing the integrand of the 2. order Cahn Hilliard problem.
*/

class CahnHilliard : public IntegrandBase
{
public:
  //! \brief The constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  CahnHilliard(unsigned short int n);
  //! \brief Empty destructor.
  virtual ~CahnHilliard() {}

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  virtual void initIntegration(size_t nGp, size_t);

  //! \brief Returns that this integrand has no explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return false; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Returns the number of primary/secondary solution field components.
  virtual size_t getNoFields(int fld) const { return fld > 1 ? 2 : 1; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual std::string getField1Name(size_t, const char* prefix) const;
  //! \brief Returns the name of the secondary solution field.
  //! \param[in] prefix Name prefix
  virtual std::string getField2Name(size_t, const char* prefix) const;

  //! \brief Sets the pointer to the tensile energy buffer.
  void setTensileEnergy(const RealArray* tens) { tensileEnergy = tens; }

  //! \brief Clears the initial crack function (after first iteration).
  void clearInitialCrack() { delete initial_crack; initial_crack = nullptr; }

protected:
  double Gc;       //!< Fracture energy density
  double smearing; //!< Smearing factor in crack
  double maxCrack; //!< Maximum value in initial crack
  double stabk;    //!< Stabilization parameter
  double scale2nd; //!< Scaling factor in front of second order term

private:
  RealFunc*         initial_crack; //!< For generating initial history field
  const RealArray*  tensileEnergy; //!< Tensile energy from elasticity solver
  mutable RealArray historyField;  //!< History field for tensile energy
};


/*!
  \brief Class representing the integrand of the 4. order Cahn Hilliard problem.
*/

class CahnHilliard4 : public CahnHilliard
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of spatial dimensions
  CahnHilliard4(unsigned short int n) : CahnHilliard(n) { scale2nd = 2.0; }
  //! \brief Empty destructor.
  virtual ~CahnHilliard4() {}

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return SECOND_DERIVATIVES; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;
};

#endif
