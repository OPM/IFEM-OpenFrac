// $Id$
//==============================================================================
//!
//! \file CHMaterial.C
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Cahn-Hilliard material models.
//!
//==============================================================================


#include "CHMaterial.h"

#include "IFEM.h"
#include "tinyxml.h"
#include "Utilities.h"


void CHMaterial::parse(const TiXmlElement* elem)
{
  utl::getAttribute(elem, "Gc", Gc);
}


void CHMaterial::printLog() const
{
  IFEM::cout << "\tConstitutive Properties: "
             << "\n\t\tCritical fracture energy density = " << Gc << std::endl;
}

double CHMaterial::getFractureEnergy(const Vec3& X) const
{
  return Gc;
}
