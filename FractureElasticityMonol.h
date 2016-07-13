// $Id$
//==============================================================================
//!
//! \file FractureElasticityMonol.h
//!
//! \date Jul 10 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementation for monolithic fracture elasticity problems.
//!
//==============================================================================

#ifndef _FRACTURE_ELASTICITY_MONOL_H
#define _FRACTURE_ELASTICITY_MONOL_H

#include "FractureElasticityVoigt.h"


/*!
  \brief Class representing the integrand of elasticity problems with fracture.
  \details This sub-class models the monolithic coupling to the phase field.
*/

class FractureElasticityMonol : public FractureElasticityVoigt
{
public:
  //! \brief The constructor invokes the parent class constructor only.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] ord Order of the phase field (2 or 4)
  FractureElasticityMonol(unsigned short int n, int ord = 2);
  //! \brief Empty destructor.
  virtual ~FractureElasticityMonol() {}

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  using FractureElasticityVoigt::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  using FractureElasticityVoigt::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  using FractureElasticityVoigt::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using FractureElasticityVoigt::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld) const;
  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix) const;

protected:
  //! \brief Extracts element solution vectors from the patch solution vectors.
  //! \param[out] eV Element solution vectors
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  bool getSolution(Vectors& eV, const std::vector<int>& MNPC) const;

private:
  unsigned short int eAcc; //!< Zero-based index to element phase-field matrix
  unsigned short int eBc;  //!< Zero-based index to element phase-field vector

  double Gc;       //!< Fracture energy density
  double smearing; //!< Smearing factor in crack
  bool   use4th;   //!< If \e true, use 4th order phase field model
  double gammaInv; //!< Penalty parameter enforcing crack irreversibility
  double crtol;    //!< Phase-field treshold for irreversibility enforcement
};

#endif
