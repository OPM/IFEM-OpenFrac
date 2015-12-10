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
#define _FRACTURE_ELASTICITY_VOIGH_H

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
  //! \brief Empty destructor.
  virtual ~FractureElasticityVoigt() {}

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

protected:
  //! \brief Evaluates the stress tensor and tensile energy at current point.
  virtual bool evalStress(double lambda, double mu, double Gc,
                          const SymmTensor& epsilon,
                          double& Phi, SymmTensor& sigma) const;

  //! \brief Evaluates the stress tensor and its derivative w.r.t. the strains.
  bool evalStress(double lambda, double mu, double Gc,
                  const SymmTensor& epsilon, double& Phi,
                  SymmTensor& sigma, Matrix* dSdE) const;
};

#endif
