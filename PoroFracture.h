// $Id$
//==============================================================================
//!
//! \file PoroFracture.h
//!
//! \date Apr 15 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for poroelasticity problems with fracture.
//!
//==============================================================================

#ifndef _PORO_FRACTURE_H
#define _PORO_FRACTURE_H

#include "PoroElasticity.h"

class FractureElasticity;


/*!
  \brief Class representing the integrand of poroelasticity with fracture.
  \details This class inherits PoroElasticity and uses elements from
  FractureElasticity in addition through a private member.
*/

class PoroFracture : public PoroElasticity
{
public:
  //! \brief The constructor allocates the internal FractureElasticy object.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] m If \e true, a mixed formulation is used
  explicit PoroFracture(unsigned short int n, bool m = false);
  //! \brief The destructor deletes the internal FractureElasticy object.
  virtual ~PoroFracture();

  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  using PoroElasticity::parseMatProp;
  //! \brief Parses material properties from an XML-element.
  virtual Material* parseMatProp(const tinyxml2::XMLElement* elem);

  //! Defines the material properties.
  virtual void setMaterial(Material* mat);

  //! \brief Initializes a time integration parameter for the integrand.
  //! \param[in] i Index of the integration parameter to define
  //! \param[in] prm The parameter value to assign
  virtual void setIntegrationPrm(unsigned short int i, double prm);

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  using PoroElasticity::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch level
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes,
                           LocalIntegral& elmInt);

  using PoroElasticity::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Returns a pointer to the Gauss-point tensile energy array.
  virtual const RealArray* getTensileEnergy() const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld Which field set to consider (1=primary,2=secondary)
  virtual size_t getNoFields(int fld) const;

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix) const;

  //! \brief Returns the applied pressure in the crack.
  RealFunc* getCrackPressure() const;

protected:
  //! \brief Computes the elasticity matrices at a quadrature point.
  //! \param elMat The element matrix object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalElasticityMatrices(ElmMats& elMat, const Matrix&,
                                      const FiniteElement& fe,
                                      const Vec3& X) const;

  //! \brief Evaluates the permeability tensor at a quadrature point.
  //! \param[out] K The permeability tensor
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool formPermeabilityTensor(SymmTensor& K,
                                      const Vectors& eV,
                                      const FiniteElement& fe,
                                      const Vec3& X) const;

  //! \brief Computes the permeability tensor of the broken material.
  //! \param[out] Kcrack Permeability tensor of the broken material
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \return Estimated crack opening
  double formCrackedPermeabilityTensor(SymmTensor& Kcrack,
                                       const Vectors& eV,
                                       const FiniteElement& fe,
                                       const Vec3& X) const;

private:
  FractureElasticity* fracEl; //!< Integrand for tangent stiffness evaluation

  double L_per; //!< Characteristic length of crack-perpendicular line
  double d_min; //!< Crack phase field value for incorporating Poiseuille flow
  double Kc;    //!< Spatial permeability in fracture
  double eps;   //!< Permeability transition exponent
};

#endif
