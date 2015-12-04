// $Id$
//==============================================================================
//!
//! \file FractureElasticity.h
//!
//! \date Oct 27 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for elasticity problems with fracture.
//!
//==============================================================================

#ifndef _FRACTURE_ELASTICITY_H
#define _FRACTURE_ELASTICITY_H

#include "Elasticity.h"

class Tensor4;


/*!
  \brief Class representing the integrand of elasticity problems with fracture.
*/

class FractureElasticity : public Elasticity
{
public:
  //! \brief The constructor invokes the parent class constructor only.
  //! \param[in] n Number of spatial dimensions
  FractureElasticity(unsigned short int n);
  //! \brief The destructor deletes the internal tensile energy buffer.
  virtual ~FractureElasticity() { delete[] myPhi; }

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  virtual void initIntegration(size_t nGp, size_t);

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns a pointer to the Gauss-point tensile energy array.
  const double* getTensileEnergy() const { return myPhi; }

private:
  //! \brief Helper method for computing the fourth-order tensor \a Gab.
  void getG(const Tensor& Ma, const Tensor& Mb, Tensor4& Gab, double eps) const;
  //! \brief Helper method for computing the fourth-order tensor \a Qa.
  void getQ(const Tensor& Ma, Tensor4& Qa, double C) const;
  //! \brief Evaluates the stress tensor and its derivative w.r.t. the strains.
  bool evalStress(double lambda, double mu, double Gc,
                  const SymmTensor& epsilon, double& Phi,
                  SymmTensor& sigma, Tensor4& dSdE) const;

protected:
  double  alpha;  //!< Relaxation factor for the crack phase field
  Vector  myCVec; //!< Crack phase field values at nodal points
  double* myPhi;  //!< Tensile energy density at integration points
};

#endif
