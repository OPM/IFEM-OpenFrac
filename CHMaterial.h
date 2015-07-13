// $Id$
//==============================================================================
//!
//! \file CHMaterial.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Cahn-Hilliard material models.
//!
//==============================================================================

#ifndef _PORO_MATERIAL_H
#define _PORO_MATERIAL_H

#include "Function.h"
#include "Vec3.h"
#include "Vec3Oper.h"

class TiXmlElement;


/*!
  \brief Class representing a material model for a Chan-Hilliard problem.
*/

class CHMaterial
{
public:
  //! \brief Empty constructor.
  CHMaterial() { Gc = 3.0; }

  //! \brief Empty destructor.
  ~CHMaterial() {}

  //! \brief Parses material parementers from an XML element.
  void parse(const TiXmlElement*);

  //! \brief Prints out material parameters to the log stream.
  void printLog() const;

  //! \brief Obtain fracture energy density.
  double getFractureEnergy(const Vec3& X) const;

protected:
  double Gc; //!< Fracture energy density.
};

#endif
