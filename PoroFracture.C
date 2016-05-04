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
#include "PoroMaterial.h"
#include "FractureElasticityVoigt.h"
#include "FiniteElement.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


PoroFracture::PoroFracture (unsigned short int n) : PoroElasticity(n)
{
  fracEl = new FractureElasticityVoigt(this,n);

  L_per = 0.01;
  d_min = 0.1;
  Kc    = 83.0;
  eps   = 50.0;
}


PoroFracture::~PoroFracture()
{
  delete fracEl;
}


bool PoroFracture::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"crack"))
    return this->PoroElasticity::parse(elem) & fracEl->parse(elem);

  IFEM::cout <<"\tCrack parameters:";
  if (utl::getAttribute(elem,"Kc",Kc))
    IFEM::cout <<" Kc = "<< Kc;
  if (utl::getAttribute(elem,"eps",eps))
    IFEM::cout <<" eps = "<< eps;
  if (utl::getAttribute(elem,"L_per",L_per))
    IFEM::cout <<" L_per = "<< L_per;
  if (utl::getAttribute(elem,"d_min",d_min))
    IFEM::cout <<" d_min = "<< d_min;
  IFEM::cout << std::endl;

  return true;
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


/*!
  This method calculates the anisotropic peremeability for the broken solid
  based on a Poiseuielle flow approximation of the fluid flow in the crack.
  See Section 5.5.2 (eqs. (104)-(109)) of Miehe and Maute (2015)
  "Phase field modeling of fracture in multi-physics problems. Part III."
*/

double PoroFracture::formCrackedPermeabilityTensor (SymmTensor& Kcrack,
                                                    const Vectors& eV,
                                                    const FiniteElement& fe,
                                                    const Vec3& X) const
{
  Vec3 gradD; // Evaluate the phase field value and gradient
  double d = fracEl->evalPhaseField(gradD,eV,fe.N,fe.dNdX);
  if (d < 0.0)
  {
    std::cerr <<" *** PoroFracture::formCrackedPermeabilityTensor(X = "<< X
              <<")\n     Invalid phase field value: "<< d << std::endl;
    return d;
  }
  else if (d < d_min)
  {
    // The crack does not affect the permeability tensor at this point
    Kcrack.zero();
    return 0.0;
  }

  double d2 = gradD.length2();
  if (d2 <= 0.0)
  {
    std::cerr <<" *** PoroFracture::formCrackedPermeabilityTensor(X = "<< X
              <<")\n     Zero phase field gradient: "<< gradD << std::endl;
    return -1.0;
  }

  Tensor F(nsd); // Calculate the deformation gradient
  if (!this->formDefGradient(eV.front(),fe.N,fe.dNdX,X.x,F))
    return -2.0;

  // Compute the inverse right Cauchy-Green tensor (C^-1)
  if (Kcrack.rightCauchyGreen(F).inverse() == 0.0)
    return -3.0;

  // Compute the symmetric tensor C^-1 - (C^-1*n0)otimes(C^-1*n0)
  // (the term in the bracket [] of eq. (108) in Miehe2015pfm3)
  Vec3 CigD = Kcrack*gradD; // C^-1*gradD
  for (unsigned short int j = 1; j <= nsd; j++)
    for (unsigned short int i = 1; i <= j; i++)
      Kcrack(i,j) -= CigD(i)*CigD(j)/d2;

  // Compute the perpendicular crack stretch
  // lambda = gradD*gradD / gradD*C^-1*gradD (see eq. (106))
  double lambda = sqrt(d2 / (gradD*CigD));
  double w = lambda*L_per - L_per; // Crack opening (see eq. (107))

  // Compute the permeability tensor, scale by d^eps*Kc*w^2*J (see eq. (108))
  Kcrack *= pow(d,eps)*Kc*w*w*F.det();
  return w;
}


bool PoroFracture::formPermeabilityTensor (SymmTensor& K,
                                           const Vectors& eV,
                                           const FiniteElement& fe,
                                           const Vec3& X) const
{
  if (this->formCrackedPermeabilityTensor(K,eV,fe,X) < 0.0)
    return false;

  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat) return false;

  Vec3 permeability = pmat->getPermeability(X);
  for (size_t i = 1; i <= K.dim(); i++)
    K(i,i) += permeability(i);

  return true;
}
