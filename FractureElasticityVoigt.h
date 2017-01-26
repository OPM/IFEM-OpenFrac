// $Id$
//==============================================================================
//!
//! \file FractureElasticityVoigt.h
//!
//! \date Dec 10 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for elasticity problems with fracture.
//!
//==============================================================================

#ifndef _FRACTURE_ELASTICITY_VOIGT_H
#define _FRACTURE_ELASTICITY_VOIGT_H

#include "FractureElasticity.h"


/*!
  \brief Class representing the integrand of elasticity problems with fracture.

  \details This sub-class uses the Voigt notation of stresses and strains,
  thereby restricted to symmetric problems.
*/

class FractureElasticityVoigt : public FractureElasticity
{
public:
  //! \brief The constructor invokes the parent class constructor only.
  //! \param[in] n Number of spatial dimensions
  FractureElasticityVoigt(unsigned short int n) : FractureElasticity(n) {}
  //! \brief Constructor for integrands with a parent integrand.
  //! \param parent The parent integrand of this one
  //! \param[in] n Number of spatial dimensions
  FractureElasticityVoigt(IntegrandBase* parent, unsigned short int n)
    : FractureElasticity(parent,n) {}
  //! \brief Empty destructor.
  virtual ~FractureElasticityVoigt() {}

  using FractureElasticity::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  using FractureElasticity::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol*) const;

protected:
  //! \brief Evaluates the stress tensor and tensile energy at current point.
  virtual bool evalStress(double lambda, double mu, double Gc,
                          const SymmTensor& epsilon,
                          double* Phi, SymmTensor& sigma) const;

  //! \brief Evaluates the stress tensor and its derivative w.r.t. the strains.
  bool evalStress(double lambda, double mu, double Gc,
                  const SymmTensor& epsilon, double* Phi,
                  SymmTensor* sigma, Matrix* dSdE,
                  bool postProc = false, bool printElm = false) const;

  friend class FractureElasticNorm;
};


/*!
  \brief Class representing the integrand of elasticity norms with fracture.
*/

class FractureElasticNorm : public ElasticityNorm
{
public:
  //! \brief The constructor invokes the parent class constructor only.
  //! \param[in] p The elasticity problem to evaluate norms for
  FractureElasticNorm(FractureElasticityVoigt& p);
  //! \brief Empty destructor.
  virtual ~FractureElasticNorm() {}

  using ElasticityNorm::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Returns the number of norm quantities.
  virtual size_t getNoFields(int group) const;
  //! \brief Returns the name of a norm quantity.
  virtual std::string getName(size_t, size_t j, const char*) const;

  static int dbgElm; //!< Element index for debug output
};

#endif
