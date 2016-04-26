// $Id$
//==============================================================================
//!
//! \file PoroFracture.C
//!
//! \date Apr 15 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for poroelasticity problems with fracture.
//!
//==============================================================================

#include "PoroFracture.h"
#include "FractureElasticityVoigt.h"
#include "Utilities.h"


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


Material* PoroFracture::parseMatProp (const TiXmlElement* elem, bool)
{
  this->PoroElasticity::parseMatProp(elem,true);
  fracEl->setMaterial(material);
  return material;
}


void PoroFracture::setMaterial (Material* mat)
{
  this->PoroElasticity::setMaterial(mat);
  fracEl->setMaterial(mat);
}


void PoroFracture::setMode (SIM::SolutionMode mode)
{
  this->PoroElasticity::setMode(mode);
  fracEl->setMode(mode);
  primsol.resize(fracEl->getNoSolutions());
}


void PoroFracture::initIntegration (size_t nGp, size_t nBp)
{
  fracEl->initIntegration(nGp,nBp);
}


LocalIntegral* PoroFracture::getLocalIntegral (size_t nen,
                                               size_t, bool neumann) const
{
  LocalIntegral* elmInt = this->PoroElasticity::getLocalIntegral(nen,0,neumann);
  fracEl->setVar(nsd+1);
  return elmInt;
}


LocalIntegral* PoroFracture::getLocalIntegral (const std::vector<size_t>& nen,
                                               size_t, bool neumann) const
{
  LocalIntegral* elmInt = this->PoroElasticity::getLocalIntegral(nen,0,neumann);
  fracEl->setVar(nsd);
  return elmInt;
}


bool PoroFracture::initElement (const std::vector<int>& MNPC,
                                LocalIntegral& elmInt)
{
  // Allocating three more vectors on the element level compared to global level
  // (1 = pressure, 2 = pressure rate, 3 = phase field)
  elmInt.vec.resize(primsol.size()+3);

  // Extract element displacement and pressure vectors
  if (!this->PoroElasticity::initElement(MNPC,elmInt))
    return false;

  // Extract element phase-field, velocity and acceleration vectors
  if (!fracEl->initElement(MNPC,elmInt))
    return false;

  if (primsol.size() < 2)
    return true; // Quasi-static simulation

  // For the standard (non-mixed) formulation, the displacement and pressure
  // variables are assumed stored interleaved in the element solution vector,
  // now we need to separate them (for velocities and accelerations)
  size_t uA = elmInt.vec.size() - 2; // Index to element acceleration vector
  size_t uV = uA - 1;                // Index to element velocity vector
  size_t pV = 2;                     // Index to element pressure rate vector
  Matrix eMat(nsd+1,MNPC.size());
  eMat = elmInt.vec[uV]; // Interleaved velocity vector
  elmInt.vec[pV] = eMat.getRow(nsd+1);  // Assign pressure rate vector
  elmInt.vec[uV] = eMat.expandRows(-1); // Assign velocity vector
  eMat.resize(nsd+1,MNPC.size());
  eMat = elmInt.vec[uA]; // Interleaved acceleration vector
  elmInt.vec[uA] = eMat.expandRows(-1); // Assign acceleration vector
  // We don't need the pressure acceleration

  return true;
}


bool PoroFracture::initElement (const std::vector<int>& MNPC,
                                const std::vector<size_t>& elem_sizes,
                                const std::vector<size_t>& basis_sizes,
                                LocalIntegral& elmInt)
{
  // Allocating three more vectors on the element level compared to global level
  // (1 = pressure, 2 = pressure rate, 3 = phase field)
  elmInt.vec.resize(primsol.size()+3);

  // Extract element displacement and pressure vectors
  if (!this->PoroElasticity::initElement(MNPC,elem_sizes,basis_sizes,elmInt))
    return false;

  // Extract element phase-field, velocity and acceleration vectors
  std::vector<int>::const_iterator uEnd = MNPC.begin() + elem_sizes.front();
  if (!fracEl->initElement(std::vector<int>(MNPC.begin(),uEnd),elmInt))
    return false;

  if (primsol.size() < 2)
    return true; // Quasi-static simulation

  // Extract the element level pressure rate vector
  std::vector<int> MNPCp(uEnd,MNPC.end());
  int ierr = utl::gather(MNPCp, 0,1, primsol[primsol.size()-2], elmInt.vec[2],
                         nsd*basis_sizes.front(), basis_sizes.front());
  if (ierr == 0) return true;

  std::cerr <<" *** PoroFracture::initElement: Detected "<< ierr
            <<" node numbers out of range."<< std::endl;
  return false;
}


const RealArray* PoroFracture::getTensileEnergy () const
{
  return fracEl->getTensileEnergy();
}


bool PoroFracture::evalElasticityMatrices (ElmMats& elMat, const Matrix&,
                                           const FiniteElement& fe,
                                           const Vec3& X) const
{
  // Evaluate tangent stiffness matrix and internal forces
  return fracEl->evalInt(elMat,fe,X);
}
