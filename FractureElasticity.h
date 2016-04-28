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
class RealFunc;


/*!
  \brief Class representing the integrand of elasticity problems with fracture.
*/

class FractureElasticity : public Elasticity
{
public:
  //! \brief The constructor invokes the parent class constructor only.
  //! \param[in] n Number of spatial dimensions
  explicit FractureElasticity(unsigned short int n);
  //! \brief Constructor for integrands with a parent integrand.
  //! \param parent The parent integrand of this one
  //! \param[in] n Number of spatial dimensions
  FractureElasticity(IntegrandBase* parent, unsigned short int n);
  //! \brief Empty destructor.
  virtual ~FractureElasticity() {}

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Sets the number of solution variables per node.
  void setVar(unsigned short int n) { npv = n; }

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return m_mode != SIM::RECOVERY; }

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  virtual void initIntegration(size_t nGp, size_t);
  //! \brief Initializes the integrand for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  virtual void initIntegration(const TimeDomain& prm, const Vector&, bool);

  using Elasticity::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  using Elasticity::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  using Elasticity::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

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

  //! \brief Evaluates the phase field and gradient at current point.
  //! \param[out] gradD Phase field gradient at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \return Phase field value at current point
  double evalPhaseField(Vec3& gradD, const Vectors& eV,
                        const Vector& N, const Matrix& dNdX) const;

  //! \brief Returns a pointer to the Gauss-point tensile energy array.
  virtual const RealArray* getTensileEnergy() const { return &myPhi; }

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] pfx Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* pfx) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note Not implemented for the tensor-based formulation.
  virtual NormBase* getNormIntegrand(AnaSol*) const { return nullptr; }

  //! \brief Returns the applied pressure in the crack.
  RealFunc* getCrackPressure() const { return crackP; }

protected:
  //! \brief Evaluates the stress tensor and tensile energy at current point.
  virtual bool evalStress(double lambda, double mu, double Gc,
                          const SymmTensor& epsilon, double* Phi,
                          SymmTensor& sigma) const;

  //! \brief Evaluates the stress tensor and its derivative w.r.t. the strains.
  bool evalStress(double lambda, double mu, double Gc,
                  const SymmTensor& epsilon, double* Phi,
                  SymmTensor& sigma, Tensor4* dSdE) const;

  //! \brief Evaluates the stress degradation function \a g(c) at current point.
  double getStressDegradation(const Vector& N, const Vectors& eV) const;

  //! \brief Evaluates Miehe's crack driving state function (eq. 56).
  double MieheCrit56(const Vec3& eps, double lambda, double mu) const;

  //! \brief Calculates integration point crack force vector contributions.
  //! \param ES Element vector to receive the force contributions
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool formCrackForce(Vector& ES, const Vectors& eV,
                      const FiniteElement& fe, const Vec3& X) const;

private:
  unsigned short int eC; //!< Zero-based index to element phase field vector

  RealFunc* crackP; //!< Applied pressure in the crack
  double    crpCut; //!< Phase-field cut-off for the applied crack pressure

  double alpha;  //!< Relaxation factor for the crack phase field
  double alpha0; //!< Initial relaxation factor
  double alpha1; //!< Relaxation factor for all steps except the first one
  Vector myCVec; //!< Crack phase field values at control (nodal) points

protected:
  double sigmaC; //!< Critical fracture tensile stress
  double zeta;   //!< Slope parameter for the driving crack force
  double tSplit; //!< No strain energy split before this time (< 0.0: always)

  mutable RealArray myPhi; //!< Tensile energy density at integration points
  Vectors&          mySol; //!< Primary solution vectors for current patch
};

#endif
