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
#include "AnaSol.h"
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "Vec3Oper.h"
#include "tinyxml.h"


CahnHilliard::CahnHilliard (unsigned short int n) : IntegrandBase(n),
  initial_crack(nullptr), flux(nullptr), tensileEnergy(nullptr), Lnorm(0)
{
  Gc = smearing = 1.0;
  maxCrack = 1.0e-3;
  scale2nd = 4.0;
  stabk = gammaInv = pthresh = 0.0;

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
  else if ((value = utl::getValue(elem,"penalty_factor")))
  {
    gammaInv = 1.0/atof(value);
    utl::getAttribute(elem,"threshold",pthresh);
  }
  else if ((value = utl::getValue(elem,"initial_crack")))
  {
    std::string type;
    utl::getAttribute(elem,"type",type);
    IFEM::cout <<"\tInitial crack function";
    initial_crack = utl::parseRealFunc(value,type);
    IFEM::cout << std::endl;
  }
  else if ((value = utl::getValue(elem,"Lnorm")))
    Lnorm = atoi(value);

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
  if (initial_crack && gammaInv <= 0.0)
    IFEM::cout <<"\n\tInitial crack specified as a function.";
  if (scale2nd == 2.0)
    IFEM::cout <<"\n\tUsing fourth-order phase field.";
  if (gammaInv > 0.0)
    IFEM::cout <<"\n\tEnforcing crack irreversibility using penalty formulation"
               <<"\n\t  gamma="<< 1.0/gammaInv
               <<" threshold="<< pthresh << std::endl;
  else
    IFEM::cout <<"\n\tEnforcing crack irreversibility using history buffer.";

  IFEM::cout << std::endl;
}


void CahnHilliard::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  primsol.resize(mode >= SIM::RHS_ONLY || gammaInv > 0.0 ? 1 : 0);
}


void CahnHilliard::initIntegration (size_t nIp, size_t)
{
  historyField.clear();
  historyField.resize(nIp,0.0);
}


bool CahnHilliard::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                            const Vec3& X) const
{
  double& H = historyField[fe.iGP];
  double GcOell = Gc/smearing;

  if (initial_crack && !(tensileEnergy && gammaInv > 0.0)) {
    // Initialize the history field using the specified initial crack function
    double dist = (*initial_crack)(X);
    if (dist < smearing)
      H = (0.25*GcOell) * (1.0/maxCrack-1.0) * (1.0-dist/smearing);
  }

  // Update history field
  if (tensileEnergy) {
    double Psi = (*tensileEnergy)[fe.iGP];
    if (gammaInv > 0.0 || Psi > H) H = Psi;
  }

  // Evaluate the current phase field, if provided
  double C = elmInt.vec.front().empty() ? 1.0 : fe.N.dot(elmInt.vec.front());

  double scale = 1.0 + 4.0*(1.0-stabk)*H/GcOell;
  if (gammaInv > 0.0 && C < pthresh)
    scale += gammaInv/GcOell;
  double s1JxW = scale*fe.detJxW;
  double s2JxW = scale2nd*smearing*smearing*fe.detJxW;

  Matrix& A = static_cast<ElmMats&>(elmInt).A.front();
  A.outer_product(fe.N, fe.N, true, s1JxW);
  A.multiply(fe.dNdX, fe.dNdX, false, true, true, s2JxW);

  static_cast<ElmMats&>(elmInt).b.front().add(fe.N,fe.detJxW);

  return true;
}


bool CahnHilliard::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                            const Vec3& X, const Vec3& normal) const
{
  double val = 0.0;

  if (flux)
    val = normal * (*flux)(X);

  static_cast<ElmMats&>(elmInt).b.front().add(fe.N,val*fe.detJxW);

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


NormBase* CahnHilliard::getNormIntegrand (AnaSol* a) const
{
  return new CahnHilliardNorm(*const_cast<CahnHilliard*>(this),Lnorm,a);
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


CahnHilliardNorm::CahnHilliardNorm (CahnHilliard& p, int Ln, const AnaSol* a)
  : NormBase(p), aSol(a)
{
  finalOp = ASM::NONE;
  Lnorm   = Ln;
}


bool CahnHilliardNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                const Vec3& X) const
{
  CahnHilliard& ch = static_cast<CahnHilliard&>(myProblem);
  ElmNorm&   pnorm = static_cast<ElmNorm&>(elmInt);

  double Gc = ch.getCriticalFracEnergy();
  double l0 = ch.getSmearingFactor();

  size_t k = 0;
  pnorm[k++] += fe.detJxW; // element volume
  for (size_t i = 0; i <= pnorm.psol.size(); i++)
  {
    const Vector& pvec = i == 0 ? elmInt.vec.front() : pnorm.psol[i-1];
    double C = pvec.dot(fe.N);
    Vector gradC;
    if (!fe.dNdX.multiply(pvec,gradC,true))
      return false;

    if (Lnorm == 1)
      pnorm[k] += fabs(C)*fe.detJxW; // L1-norm, |c|
    else if (Lnorm == 2)
      pnorm[k] += C*C*fe.detJxW; // L2-norm, |c|
    else if (Lnorm == -1 && C > 0.0)
      if (pnorm[k] == 0.0 || C < pnorm[k])
        pnorm[k] = C; // Smallest-value norm

    if (Lnorm > 0)
      k += 2; // Make space for the volume-specific norm |c|/V
    else if (Lnorm < 0)
      k ++;

    // Dissipated energy, eps_d
    pnorm[k++] += Gc*(pow(C-1.0,2.0)/(4.0*l0) + l0*gradC.dot(gradC))*fe.detJxW;

    if (aSol && i == 0) { // add final group
      size_t nNorms = pnorm.size();
      pnorm[nNorms-6] += C*C*fe.detJxW;
      pnorm[nNorms-5] += gradC.dot(gradC)*fe.detJxW;
      double CA = (*aSol->getScalarSol())(X);
      Vec3 gradCA = (*aSol->getScalarSecSol())(X);
      pnorm[nNorms-4] += CA*CA*fe.detJxW;
      pnorm[nNorms-3] += gradCA*gradCA*fe.detJxW;
      double eC = C-CA;
      pnorm[nNorms-2] += eC*eC*fe.detJxW;
      gradCA -= gradC;
      pnorm[nNorms-1] += gradCA*gradCA*fe.detJxW;
    }
  }

  return true;
}


bool CahnHilliardNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (Lnorm < 1) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the volume-specific norm |c|/V
  for (size_t ip = 1; ip < pnorm.size()-(aSol?6:0); ip += 3)
    pnorm[ip+1] = pnorm[ip] / pnorm[0];

  return true;
}


std::string CahnHilliardNorm::getName (size_t i, size_t j,
                                       const char* prefix) const
{
  std::string name;
  if (i == 2) {
    static const char* names[] = {"||c^h||_L2", "||c^h||_H1",
                                  "||c||_L2", "||c||_H1",
                                  "||e^h||_L2", "||e^h||_H1"};
    name = names[j-1];
  }
  else {
    if (i == 1 && j == 1)
      return "volume";
    else if (i == 1 && j > 1)
      j --;
    if (Lnorm == 0)
      j += 2;
    else if (Lnorm < 0 && j > 1)
      j ++;

    if (j == 1)
      name = "|c|";
    else if (j == 2)
      name = "|c|/V";
    else if (j == 3)
      name = "eps_d";
    else
      return this->NormBase::getName(i,j,prefix);
  }

  if (!prefix)
    return name;

  return std::string(prefix) + " " + name;
}


size_t CahnHilliardNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields() + (aSol ? 1 : 0);
  else if (group == (int)(1+this->NormBase::getNoFields()))
    return 6;
  else if (Lnorm == 0)
    return 2;
  else if (group == 1 && Lnorm > 0)
    return 4;
  else
    return 3;
}
