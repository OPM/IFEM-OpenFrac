// $Id$
//==============================================================================
//!
//! \file PoroFracture.h
//!
//! \date Apr 15 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for elasticity problems with fracture.
//!
//==============================================================================

#ifndef _PORO_FRACTURE_H
#define _PORO_FRACTURE_H

#include "PoroElasticity.h"

class FractureElasticity;


/*!
  \brief Class representing the integrand of poroelasticity with fracture.

  \details This class inherits PoroElasticity and uses elements from
  FractureElasticity through a private member.
*/

class PoroFracture : public PoroElasticity
{
public:
  //! \brief The constructor allocates the internal FractureElasticy object.
  //! \param[in] n Number of spatial dimensions
  PoroFracture(unsigned short int n);
  //! \brief The destructor deletes the internal FractureElasticy object.
  virtual ~PoroFracture();

  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  //! \brief Returns a pointer to the Gauss-point tensile energy array.
  virtual const RealArray* getTensileEnergy() const;

protected:
  //! \brief Computes the elasticity matrices for a quadrature point.
  //! \param elmInt The element matrix object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalElasticityMatrices(ElmMats& elMat, const Matrix&,
                                      const FiniteElement& fe,
                                      const Vec3& X) const;

private:
  FractureElasticity* fracEl; //!< Integrand for tangent stiffness evaluation
};

#endif
