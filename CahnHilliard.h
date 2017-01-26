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
  \brief Class representing the integrand of a 2nd order Cahn Hilliard problem.
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

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  using IntegrandBase::initIntegration;
  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  virtual void initIntegration(size_t nGp, size_t);

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using IntegrandBase::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  using IntegrandBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns the number of primary/secondary solution field components.
  virtual size_t getNoFields(int fld) const { return fld > 1 ? 2 : 1; }
  //! \brief Returns the name of the primary solution field component.
  //! \param[in] prefix Name prefix
  virtual std::string getField1Name(size_t, const char* prefix) const;
  //! \brief Returns the name of the secondary solution field component.
  //! \param[in] idx Field component index
  //! \param[in] prefix Name prefix
  virtual std::string getField2Name(size_t idx, const char* prefix) const;

  //! \brief Sets the pointer to the tensile energy buffer.
  void setTensileEnergy(const RealArray* tens) { tensileEnergy = tens; }
  //! \brief Sets the flux function associated with the Neumann boundary.
  void setFlux(VecFunc* f) { flux = f; }

  //! \brief Returns the initial crack function.
  RealFunc* initCrack() { return initial_crack; }
  //! \brief Clears the initial crack function (used after first time step).
  void clearInitialCrack() { delete initial_crack; initial_crack = nullptr; }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* a) const;

  //! \brief Returns the critical fracture energy.
  double getCriticalFracEnergy() const { return Gc; }
  //! \brief Returns the smearing factor.
  double getSmearingFactor() const { return smearing; }
  //! \brief Scale the smearing factor, for use during initial refinement cycle.
  double scaleSmearing(double s) { return smearing *= s; }

protected:
  double Gc;       //!< Fracture energy density
  double smearing; //!< Smearing factor in crack
  double maxCrack; //!< Maximum value in initial crack
  double stabk;    //!< Stabilization parameter
  double scale2nd; //!< Scaling factor in front of second order term
  double gammaInv; //!< Penalty factor (if positive) for crack irreversibility
  double pthresh;  //!< Penalty formulation phase field threshold

private:
  RealFunc*        initial_crack; //!< For generating initial history field
  VecFunc*         flux;          //!< Boundary flux function
  const RealArray* tensileEnergy; //!< Tensile energy from elasticity solver
  int              Lnorm;         //!< Which L-norm to integrate

public:
  mutable RealArray historyField; //!< History field for tensile energy
};


/*!
  \brief Class representing the integrand of a 4th order Cahn Hilliard problem.
*/

class CahnHilliard4 : public CahnHilliard
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  CahnHilliard4(unsigned short int n) : CahnHilliard(n) { scale2nd = 2.0; }
  //! \brief Empty destructor.
  virtual ~CahnHilliard4() {}

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return SECOND_DERIVATIVES; }

  using CahnHilliard::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;
};


/*!
  \brief Class representing the norms for a Cahn Hilliard problem.
*/

class CahnHilliardNorm : public NormBase
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  CahnHilliardNorm(CahnHilliard& p, int Ln, const AnaSol* a = nullptr);
  //! \brief Empty destructor.
  virtual ~CahnHilliardNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using NormBase::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //!
  //! \details This method is used to compute volume-normalized norms.
  virtual bool finalizeElement(LocalIntegral& elmInt);

  //! \brief Returns the number of norm groups or size of a specified group.
  virtual size_t getNoFields(int group) const;
  //! \brief Returns the name of a norm quantity.
  virtual std::string getName(size_t i, size_t j, const char* prefix) const;

private:
  int Lnorm; //!< Which L-norm to integrate (0: none, 1: L1-norm, 2: L2-norm)
  const AnaSol* aSol; //!< Analytical solution
};

#endif
