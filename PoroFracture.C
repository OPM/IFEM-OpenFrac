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


/*!
  \brief Class representing the integrand of elasticity problems with fracture.

  \details This sub-class overrides the getElementSolution method, to account
  for that the primary solution vector, as extracted from the patch level,
  also contains the pressure variables in addition to the displacements.
*/

class PoroFractureElasticity : public FractureElasticityVoigt
{
public:
  //! \brief Constructor for integrands with a parent integrand.
  //! \param parent The parent integrand of this one
  //! \param[in] n Number of spatial dimensions
  PoroFractureElasticity(IntegrandBase* parent, unsigned short int n)
    : FractureElasticityVoigt(parent,n) {}
  //! \brief Empty destructor.
  virtual ~PoroFractureElasticity() {}

  //! \brief Retrieves the element solution vectors.
  //! \param[out] eV Element solution vectors
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool getElementSolution(Vectors& eV,
                                  const std::vector<int>& MNPC) const
  {
    if (!this->FractureElasticityVoigt::getElementSolution(eV,MNPC))
      return false;

    // Filter out the pressure components
    // FIXME: Mixed
    size_t nen = MNPC.size();
    if (eV.front().size() == nen*(nsd+1))
    {
      Vector temp(eV.front());
      Vector& actual = eV.front();
      actual.resize(nen*nsd);
      for (size_t n = 0; n < nen; n++)
        for (unsigned short int i = 0; i < nsd; i++)
          actual[nsd*n+i] = temp[(nsd+1)*n+i];
    }

    return true;
  }
};


PoroFracture::PoroFracture (unsigned short int n) : PoroElasticity(n)
{
  fracEl = new PoroFractureElasticity(this,n);

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

  // Extract element displacement, velocity, acceleration and pressure vectors
  if (!this->PoroElasticity::initElement(MNPC,elmInt))
    return false;

  // Extract element phase-field vector
  return fracEl->initElement(MNPC,elmInt);
}


bool PoroFracture::initElement (const std::vector<int>& MNPC,
                                const std::vector<size_t>& elem_sizes,
                                const std::vector<size_t>& basis_sizes,
                                LocalIntegral& elmInt)
{
  // Allocating three more vectors on the element level compared to global level
  // (1 = pressure, 2 = pressure rate, 3 = phase field)
  elmInt.vec.resize(primsol.size()+3);

  // Extract element displacement, velocity, acceleration and pressure vectors
  if (!this->PoroElasticity::initElement(MNPC,elem_sizes,basis_sizes,elmInt))
    return false;

  // Extract element phase-field vector
  std::vector<int>::const_iterator uEnd = MNPC.begin() + elem_sizes.front();
  return fracEl->initElement(std::vector<int>(MNPC.begin(),uEnd),elmInt);
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
  This method calculates the anisotropic permeability for the broken solid
  based on a Poiseuille flow approximation of the fluid flow in the crack.
  See Section 5.5.2 (eqs. (104)-(109)) of Miehe and Maute (2015)
  "Phase field modeling of fracture in multi-physics problems. Part III."
  as well as clarifying note by Fonn and Kvamsdal (2016).
*/

double PoroFracture::formCrackedPermeabilityTensor (SymmTensor& Kcrack,
                                                    const Vectors& eV,
                                                    const FiniteElement& fe,
                                                    const Vec3& X) const
{
  // Base permeability
  // FIXME: What to do for non-isotropic permeability?
  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat)
    return -4.0;
  double perm = pmat->getPermeability(X)[0];

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
    Kcrack = perm;
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

  // Compute alpha = |grad D| / |F^-T grad D|, see F&K eq. (44)
  Tensor Finv(F, true);
  Finv.inverse();
  Vec3 FtgradD = Finv * gradD;
  double alpha = gradD.length() / FtgradD.length();
  SymmTensor kCinv = perm * Kcrack; // k * C^-1

  // Compute the symmetric tensor C^-1 - alpha * (C^-1*n0)otimes(C^-1*n0)
  // (the term in the bracket [] of eq. (108) in Miehe2015pfm3)
  // See also F&K eq. (45)
  Vec3 CigD = Kcrack*gradD; // C^-1*gradD
  double a2od2 = alpha*alpha/d2;
  for (unsigned short int j = 1; j <= nsd; j++)
    for (unsigned short int i = 1; i <= j; i++)
      Kcrack(i,j) -= a2od2*CigD(i)*CigD(j);

  // Compute the perpendicular crack stretch
  // lambda = gradD*gradD / gradD*C^-1*gradD (see M&M eq. (106), F&K eq. (36))
  double lambda = sqrt(d2 / (gradD*CigD));
#if INT_DEBUG > 3
  std::cout <<"PoroFracture::formCrackedPermeabilityTensor(X = "<< X
            <<") lambda = "<< lambda << std::endl;
#endif
  if (lambda <= 1.0)
  {
    Kcrack = perm;
    return 0.0;
  }

  // Compute the permeability tensor, scale by d^eps*Kc*w^2*J (see eq. (108))
  // See also F&K eq. (45)
  double w = lambda*L_per - L_per; // Crack opening (see M&M eq. (107))
  if (w < 0.0) w = 0.0;            // See F&K eq. (37)
  Kcrack *= pow(d,eps) * (w*w/12.0 - perm);
  Kcrack += kCinv;
  Kcrack *= F.det();

  return w;
}


bool PoroFracture::formPermeabilityTensor (SymmTensor& K,
                                           const Vectors& eV,
                                           const FiniteElement& fe,
                                           const Vec3& X) const
{
  return this->formCrackedPermeabilityTensor(K,eV,fe,X) >= 0.0;
}


size_t PoroFracture::getNoFields (int fld) const
{
  size_t nK = fld == 2 ? 1+nsd : 0;
  return this->PoroElasticity::getNoFields(fld) + nK;
}


std::string PoroFracture::getField2Name (size_t i, const char* prefix) const
{
  if (i < this->PoroElasticity::getNoFields(2))
    return this->PoroElasticity::getField2Name(i,prefix);

  i -= this->PoroElasticity::getNoFields(2);

  static const char* s[4] = { "w", "K_xx", "K_yy", "K_zz" };

  if (!prefix) return s[i%4];

  return prefix + std::string(" ") + s[i%4];
}


bool PoroFracture::evalSol (Vector& s, const FiniteElement& fe,
                            const Vec3& X, const std::vector<int>& MNPC) const
{
  if (!this->PoroElasticity::evalSol(s,fe,X,MNPC))
    return false;

  Vectors eV;
  if (!fracEl->getElementSolution(eV,MNPC))
    return false;

  SymmTensor Kc(nsd);
  s.push_back(this->formCrackedPermeabilityTensor(Kc,eV,fe,X));
  for (size_t i = 1; i <= Kc.dim(); i++)
    s.push_back(Kc(i,i));

  return true;
}
