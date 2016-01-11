// $Id$
//==============================================================================
//!
//! \file CahnHillard.C
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Cahn-Hilliard problems.
//!
//==============================================================================

#include "CahnHilliard.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


CahnHilliard::CahnHilliard (unsigned short int n) : IntegrandBase(n),
  Gc(1.0), smearing(1.0), maxCrack(1.0e-3), stabk(0.0), scale2nd(4.0),
  initial_crack(nullptr), tensileEnergy(nullptr)
{
  primsol.resize(1);
}


bool CahnHilliard::parse (const TiXmlElement* elem)
{
  const char* value = utl::getValue(elem,"Gc");
  if (value)
    Gc = atof(value);
  else if ((value = utl::getValue(elem,"smearing")))
    smearing = atof(value);
  else if ((value = utl::getValue(elem,"maxcrack")))
    maxCrack = atof(value);
  else if ((value = utl::getValue(elem,"stabilization")))
    stabk = atof(value);
  else if ((value = utl::getValue(elem,"initial_crack")))
  {
    std::string type;
    utl::getAttribute(elem,"type",type);
    IFEM::cout <<"\tInitial crack function";
    initial_crack = utl::parseRealFunc(value,type);
    IFEM::cout << std::endl;
  }

  return true;
}


void CahnHilliard::printLog () const
{
  IFEM::cout <<"Cahn-Hilliard: "<< nsd <<"D"
             <<"\n\tCritical fracture energy density: "<< Gc
             <<"\n\tSmearing factor: "<< smearing
             <<"\n\tMax value in crack: "<< maxCrack;
  if (stabk != 0.0)
    IFEM::cout <<"\n\tStabilization parameter: "<< stabk;
  if (initial_crack)
    IFEM::cout <<"\n\tInitial crack specified as a function.";
  if (scale2nd == 2.0)
    IFEM::cout <<"\n\tUsing fourth-order phase field.";
  IFEM::cout << std::endl;
}


void CahnHilliard::initIntegration (size_t nIp, size_t)
{
  historyField.resize(nIp,0.0);
}


bool CahnHilliard::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                            const Vec3& X) const
{
  double& H = historyField[fe.iGP];

  if (initial_crack) {
    double dist = (*initial_crack)(X);
    if (dist < smearing)
      H = (0.25*Gc/smearing) * (1.0/maxCrack-1.0) * (1.0-dist/smearing);
  }

  // Update history field
  if (tensileEnergy)
    H = std::max(H,(*tensileEnergy)[fe.iGP]);

  double scale = 1.0 + 4.0*smearing*(1.0-stabk)*H/Gc;
  double s1JxW = scale*fe.detJxW;
  double s2JxW = scale2nd*smearing*smearing*fe.detJxW;

  Matrix& A = static_cast<ElmMats&>(elmInt).A.front();
  for (size_t i = 1; i <= fe.N.size(); i++)
    for (size_t j = 1; j <= fe.N.size(); j++) {
      double grad = 0.0;
      for (size_t k = 1; k <= nsd; k++)
        grad += fe.dNdX(i,k)*fe.dNdX(j,k);
      A(i,j) += fe.N(i)*fe.N(j)*s1JxW + grad*s2JxW;
    }

  static_cast<ElmMats&>(elmInt).b.front().add(fe.N,fe.detJxW);

  return true;
}


bool CahnHilliard::evalSol (Vector& s, const FiniteElement& fe,
                            const Vec3& X, const std::vector<int>& MNPC) const
{
  s.resize(2);

  Vector tmp;
  if (utl::gather(MNPC,1,primsol.front(),tmp))
    return false;

  double c = fe.N.dot(tmp);
  if (c < maxCrack)
    s(1) = 0.0;
  else if (c > 1.0)
    s(1) = 1.0;
  else
    s(1) = c;

  s(2) = historyField[fe.iGP];
  return true;
}


std::string CahnHilliard::getField1Name (size_t, const char* prefix) const
{
  return prefix ? prefix + std::string(" phase") : std::string("phase");
}


std::string CahnHilliard::getField2Name (size_t idx, const char* prefix) const
{
  std::string name(idx == 0 ? "phase" : "Max tensile energy");
  return prefix ? std::string(prefix) + " " + name : name;
}


NormBase* CahnHilliard::getNormIntegrand (AnaSol*) const
{
  return new CahnHilliardNorm(*const_cast<CahnHilliard*>(this));
}


bool CahnHilliard4::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                             const Vec3& X) const
{
  if (!this->CahnHilliard::evalInt(elmInt,fe,X))
    return false;

  Matrix& A = static_cast<ElmMats&>(elmInt).A.front();
  double s4JxW = pow(smearing,4.0)*fe.detJxW;

  for (size_t i = 1; i <= fe.N.size(); i++)
    for (size_t j = 1; j <= fe.N.size(); j++) {
      double grad = 0.0;
      for (unsigned short int k = 1; k <= nsd; k++)
        grad += fe.d2NdX2(i,k,k)*fe.d2NdX2(j,k,k);
      A(i,j) += grad*s4JxW;
    }

  return true;
}


bool CahnHilliardNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                const Vec3& X) const
{
  CahnHilliard& ch = static_cast<CahnHilliard&>(myProblem);
  ElmNorm&   pnorm = static_cast<ElmNorm&>(elmInt);

  double Gc = ch.getCriticalFracEnergy();
  double l0 = ch.getSmearingFactor();

  size_t k = 0;
  for (size_t i = 0; i <= pnorm.psol.size(); i++) {
    const Vector& pvec = i == 0 ? elmInt.vec.front() : pnorm.psol[i-1];
    double C = pvec.dot(fe.N);
    Vector gradC;
    if (!fe.dNdX.multiply(pvec,gradC,true))
      return false;

    pnorm[k++] += C*C*fe.detJxW;
    pnorm[k++] += Gc*(pow(C-1.0,2.0)/(4.0*l0) + l0*gradC.dot(gradC))*fe.detJxW;
  }

  return true;
}


size_t CahnHilliardNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else
    return 2;
}
