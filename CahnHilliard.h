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

#ifndef _CAHNHILLIARD_H
#define _CAHNHILLIARD_H

#include "IntegrandBase.h"
#include "Vec3.h"
#include "CHMaterial.h"


/*!
  \brief Class representing the integrand of the 2. order Cahn Hilliard problem.
*/

class CahnHilliard : public IntegrandBase
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  CahnHilliard(unsigned short int n = 3);
  //! \brief The destructor deletes the functions to be Galerkin-projected.
  virtual ~CahnHilliard() {}

  //! \brief Set smearing factor.
  void setSmearFactor(double smear) { smearFactor = smear; }

  //! \brief Set cap value in crack.
  void setMaxCrack(double crack) { maxCrack = crack; }

  //! \brief Returns current material.
  const CHMaterial* getMaterial() const { return mat; }

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

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
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const { return fld > 1 ? 1 : 1; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual std::string getField1Name(size_t, const char* prefix = nullptr) const;

  //! \brief Set material parameters for current patch.
  void setMaterial(CHMaterial* material) { mat = material; }

  //! \brief Set a function for the initial crack. Used to generate history field.
  void setInitialCrackFunction(RealFunc* func) { initial_crack = func; }

  //! \brief Set tensile energy vector.
  void setTensileEnergy(const Vector* tens) { tensile = tens; }

  //!< Print problem definition to log.
  virtual void printLog();

protected:
  unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)
  double smearFactor; //!< Smearing factor in crack.
  double maxCrack; //!< Maximum value in initial crack.
  double stabk; //!< Stabilization parameter.
  // Need to do this here as we do not have access to integration point coordinates on SIM level.
  RealFunc*  initial_crack; //!< For generating the initial history field on first integration.

  CHMaterial* mat; //!< Material parameters.

  mutable Vector historyField;   //!< History field for tensile energy.
  const Vector* tensile; //!< Tensile energy from elasticity solver.

  double second_scale; //!< Scaling factor in front of second order term
};


/*!
  \brief Class representing the integrand of the 4. order Cahn Hilliard problem.
*/

class CahnHilliard4 : public CahnHilliard
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of spatial dimensions
  CahnHilliard4(unsigned short int n = 3);

  //! \brief Empty destructor.
  virtual ~CahnHilliard4() {}

  //! \brief Returns integrand traits.
  virtual int getIntegrandType() const { return SECOND_DERIVATIVES; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;
};

#endif
