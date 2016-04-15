// $Id$
//==============================================================================
//!
//! \file PoroFracture.C
//!
//! \date Apr 15 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for elasticity problems with fracture.
//!
//==============================================================================

#include "PoroFracture.h"
#include "FractureElasticityVoigt.h"


PoroFracture::PoroFracture (unsigned short int n) : PoroElasticity(n)
{
  fracEl = new FractureElasticityVoigt(this,n);
}


PoroFracture::~PoroFracture()
{
  delete fracEl;
}


bool PoroFracture::parse (const TiXmlElement* elem)
{
  return this->PoroElasticity::parse(elem) & fracEl->parse(elem);
}


void PoroFracture::initIntegration (size_t nGp, size_t nBp)
{
  fracEl->initIntegration(nGp,nBp);
}


bool PoroFracture::initElement (const std::vector<int>& MNPC,
                                LocalIntegral& elmInt)
{
  return this->PoroElasticity::initElement(MNPC,elmInt) &&
         fracEl->initElement(MNPC,elmInt);
}


const RealArray* PoroFracture::getTensileEnergy () const
{
  return fracEl->getTensileEnergy();
}


bool PoroFracture::evalElasticityMatrices (ElmMats& elMat, const Matrix&,
                                           const FiniteElement& fe,
                                           const Vec3& X) const
{
  return fracEl->evalInt(elMat,fe,X);
}
