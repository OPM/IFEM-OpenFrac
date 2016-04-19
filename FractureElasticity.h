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
  //! \brief Constructor for integrands with a parent integrand.
  //! \param parent The parent integrand of this one
  //! \param[in] n Number of spatial dimensions
  FractureElasticity(IntegrandBase* parent, unsigned short int n);
  //! \brief Empty destructor.
  virtual ~FractureElasticity() {}

  //! \brief Sets the number of solution variables per node.
  void setVar(unsigned short int n) { npv = n; }

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

  using Elasticity::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  //! \param[out] pdir Directions of the principal stresses
  virtual bool evalSol(Vector& s, const Vectors& eV, const FiniteElement& fe,
                       const Vec3& X, bool toLocal, Vec3* pdir) const;

  //! \brief Returns a pointer to the Gauss-point tensile energy array.
  const RealArray* getTensileEnergy() const { return &myPhi; }

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* pfx) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note Not implemented for the tensor-based formulation.
  virtual NormBase* getNormIntegrand(AnaSol*) const { return nullptr; }

protected:
  //! \brief Evaluates the stress tensor and tensile energy at current point.
  virtual bool evalStress(double lambda, double mu, double Gc,
                          const SymmTensor& epsilon,
                          double* Phi, SymmTensor& sigma) const;

  //! \brief Evaluates the stress tensor and its derivative w.r.t. the strains.
  bool evalStress(double lambda, double mu, double Gc,
                  const SymmTensor& epsilon, double* Phi,
                  SymmTensor& sigma, Tensor4* dSdE,
                  bool postProc = false) const;

  //! \brief Evaluates the stress degradation function \a g(c) at current point.
  double getStressDegradation(const Vector& N, const Vectors& eV) const;

private:
  unsigned short int eC; //!< Zero-based index to element phase field vector

protected:
  double alpha;  //!< Relaxation factor for the crack phase field
  Vector myCVec; //!< Crack phase field values at nodal points

  mutable RealArray myPhi; //!< Tensile energy density at integration points
  Vectors&          mySol; //!< Primary solution vectors for current patch
};

#endif
