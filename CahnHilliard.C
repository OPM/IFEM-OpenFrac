// $Id$
//==============================================================================
//!
//! \file CahnHilliard.C
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
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"


CahnHilliard::CahnHilliard (unsigned short int n) : IntegrandBase(n),
  initial_crack(nullptr), flux(nullptr), tensileEnergy(nullptr), Lnorm(0)
{
  Gc = smearing = l0 = 1.0;
  maxCrack = 1.0e-3;
  scale2nd = 4.0;
  stabk = gammaInv = pthresh = 0.0;
}


bool CahnHilliard::parse (const TiXmlElement* elem)
{
  const char* value = utl::getValue(elem,"Gc");
  if (value)
    Gc = atof(value);
  else if ((value = utl::getValue(elem,"smearing")))
  {
    smearing = atof(value);
    utl::getAttribute(elem,"l0",l0);
    if (l0 > smearing) l0 = smearing;
  }
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
  if (initial_crack && gammaInv == 0.0)
    IFEM::cout <<"\n\tInitial crack specified as a function.";
  if (scale2nd == 2.0)
    IFEM::cout <<"\n\tUsing fourth-order phase field.";
  if (gammaInv != 0.0)
    IFEM::cout <<"\n\tEnforcing crack irreversibility using penalty formulation"
               <<"\n\t  gamma="<< 1.0/(gammaInv > 0.0 ? gammaInv : -gammaInv)
               <<" threshold="<< pthresh << std::endl;
  else
    IFEM::cout <<"\n\tEnforcing crack irreversibility using history buffer.";

  IFEM::cout << std::endl;
}


void CahnHilliard::clearInitialCrack ()
{
  delete initial_crack;
  initial_crack = nullptr;
}


void CahnHilliard::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  primsol.resize(gammaInv != 0.0 ? 2 : (mode < SIM::RHS_ONLY ? 0 : 1));
}


void CahnHilliard::initIntegration (size_t nIp, size_t)
{
  historyField.clear();
  historyField.resize(nIp,0.0);
}


LocalIntegral* CahnHilliard::getLocalIntegral (size_t nen, size_t,
                                               bool neumann) const
{
  LocalIntegral* li = this->IntegrandBase::getLocalIntegral(nen,0,neumann);
  if (m_mode >= SIM::RHS_ONLY && !neumann)
    static_cast<ElmMats*>(li)->c.resize(1); // Total dissipation energy
  return li;
}


bool CahnHilliard::evalIntD (ElmMats& elm, const FiniteElement& fe) const
{
  // Store the tensile energy density in the history buffer
  historyField[fe.iGP] = (*tensileEnergy)[fe.iGP];

  // Evaluate the previous phase field, if provided
  double C = elm.vec.back().empty() ? 1.0 : fe.N.dot(elm.vec.back());
  bool inCrack = C < pthresh;
#if INT_DEBUG > 3
  std::cout <<"\nCahnHilliard::evalIntD("<< fe.iGP <<"): C = "<< C;
#endif

  double GcOell = 0.5*Gc/smearing; // Note: ell = 2*smearing
  double scale = GcOell + 2.0*(1.0-stabk)*historyField[fe.iGP];
  if (inCrack) scale -= gammaInv; // Note: gammaInv is assumed negative here
  double s1JxW = scale*fe.detJxW;
  double s2JxW = scale2nd*0.5*Gc*smearing*fe.detJxW;
#if INT_DEBUG > 3
  if (inCrack)
    std::cout <<"\n\tIn crack: scale "<< scale+gammaInv <<" --> "<< scale;
  std::cout << std::endl;
#endif

  if (m_mode == SIM::STATIC)
  {
    Matrix& A = elm.A.front();

    A.outer_product(fe.N,fe.N,true,s1JxW);             // A +=  N  * N^t  *s1JxW
    A.multiply(fe.dNdX,fe.dNdX,false,true,true,s2JxW); // A += dNdX*dNdX^t*s2JxW

    double s3JxW = (scale - GcOell)*fe.detJxW;
    elm.b.front().add(fe.N,s3JxW); // R += N*s3JxW
  }
  else if (m_mode == SIM::INT_FORCES && !elm.vec.front().empty())
  {
    Vector& R = elm.b.front();
    double& E = elm.c.front();

    // Evaluate the current phase field.
    // Note that we do not cap the value to fit within the range [0,1] here,
    // because it needs to be consistent with the solution itself.
    C = fe.N.dot(elm.vec.front());

    Vector gradC; // Evaluate the phase field gradient gradC = dNdX^t*eC
    if (!fe.dNdX.multiply(elm.vec.front(),gradC,true))
      return false;

    // Integrate the dissipated energy.
    // Note that the penalty term is also included here,
    // in contrast to for the CahnHilliardNorm integrand.
    // Therefore the values will appear different.
    E += (0.5*GcOell*(1.0-C)*(1.0-C) + Gc*smearing*gradC.dot(gradC))*fe.detJxW;
    if (inCrack)
      E -= 0.5*gammaInv*C*C*fe.detJxW;

#if INT_DEBUG > 3
    std::cout <<"\tC = "<< C <<"  E = "<< E << std::endl;
#endif

    // Integrate the residual force vector
    double s3JxW = (GcOell - scale*C)*fe.detJxW;

    R.add(fe.N,s3JxW); // R += N*s3JxW

    gradC *= s2JxW;
    return fe.dNdX.multiply(gradC,R,false,true); // R += dNdX*gradC*s2JxW
  }
  else
  {
    std::cerr <<" *** CahnHilliard::evalIntD: Invalid simulation mode "
              << m_mode << std::endl;
    return false;
  }

  return true;
}


bool CahnHilliard::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                            const Vec3& X) const
{
  if (this->useDformulation()) // d=1-c is to be the primary unknown
    return this->evalIntD(static_cast<ElmMats&>(elmInt),fe);

  double& H = historyField[fe.iGP];
  double GcOell = 0.5*Gc/smearing; // Note: ell = 2*smearing
  double dist = 0.0, Psi = 0.0;

  // Initialize the history field using the specified initial crack function
  if (initial_crack && !(tensileEnergy && gammaInv > 0.0))
    if ((dist = (*initial_crack)(X)) < smearing)
      H = (0.5*GcOell) * (1.0/maxCrack-1.0) * (1.0-dist/smearing);

  // Update history field
  if (tensileEnergy && tensileEnergy->size() == historyField.size())
    if ((Psi = (*tensileEnergy)[fe.iGP]) > H || gammaInv > 0.0)
      H = Psi;

  // Evaluate the previous phase field, if provided
  double C = elmInt.vec.back().empty() ? 1.0 : fe.N.dot(elmInt.vec.back());
  bool inCrack = gammaInv > 0.0 && C < pthresh;
#if INT_DEBUG > 3
  std::cout <<"\nCahnHilliard::evalInt(X = "<< X <<"): C = "<< C;
  if (Psi > 0.0 || (dist > 0.0 && dist < smearing))
  {
    std::cout <<"\n\tHistory field, iGp="<< fe.iGP;
    if (Psi > 0.0)
      std::cout <<": (from tensile energy = "<< Psi;
    else
      std::cout <<": (from initial crack = "<< dist;
    std::cout <<") "<< H;
  }
#endif

  double scale = 1.0 + 2.0*(1.0-stabk)*H/GcOell;
  if (inCrack)
    scale += gammaInv/GcOell;
  double s1JxW = scale*fe.detJxW;
  double s2JxW = scale2nd*smearing*smearing*fe.detJxW;
#if INT_DEBUG > 3
  if (inCrack)
    std::cout <<"\n\tIn crack: scale "<< scale-gammaInv/GcOell <<" -> "<< scale;
  std::cout << std::endl;
#endif

  if (m_mode == SIM::STATIC)
  {
    Matrix& A = static_cast<ElmMats&>(elmInt).A.front();

    A.outer_product(fe.N,fe.N,true,s1JxW);             // A +=  N  * N^t  *s1JxW
    A.multiply(fe.dNdX,fe.dNdX,false,true,true,s2JxW); // A += dNdX*dNdX^t*s2JxW

    static_cast<ElmMats&>(elmInt).b.front().add(fe.N,fe.detJxW); // R += N*|J|*W
  }
  else if (m_mode == SIM::INT_FORCES && !elmInt.vec.front().empty())
  {
    Vector& R = static_cast<ElmMats&>(elmInt).b.front();
    double& E = static_cast<ElmMats&>(elmInt).c.front();

    // Evaluate the current phase field.
    // Note that we do not cap the value to fit within the range [0,1] here,
    // because it needs to be consistent with the solution itself.
    C = fe.N.dot(elmInt.vec.front());

    Vector gradC; // Compute the phase field gradient gradC = dNdX^t*eC
    if (!fe.dNdX.multiply(elmInt.vec.front(),gradC,true))
      return false;

    // Integrate the dissipated energy.
    // Note that the penalty term is also included here,
    // in contrast to for the CahnHilliardNorm integrand.
    // Therefore the values will appear different.
    E += (0.5*GcOell*(1.0-C)*(1.0-C) + Gc*smearing*gradC.dot(gradC))*fe.detJxW;
    if (inCrack)
      E += 0.5*gammaInv*C*C*fe.detJxW;

#if INT_DEBUG > 3
    std::cout <<"\tC = "<< C <<"  E = "<< E << std::endl;
#endif

    // Integrate the residual force vector.
    // Apply scaling Gc/ell compared to the STATIC mode, such that
    // the resulting residual force vector has comparable dimension
    // as the residual forces of the elasticity equation.
    s1JxW  = GcOell*(1.0 - C*scale)*fe.detJxW;
    s2JxW *= GcOell;

    R.add(fe.N,s1JxW); // R += N*s1JxW

    gradC *= -s2JxW;
    return fe.dNdX.multiply(gradC,R,false,true); // R -= dNdX*gradC*s2JxW
  }
  else
  {
    std::cerr <<" *** CahnHilliard::evalInt: Invalid simulation mode "
              << m_mode << std::endl;
    return false;
  }

  return true;
}


bool CahnHilliard::evalIntMx (LocalIntegral& elmInt, const MxFiniteElement& fe,
                              const Vec3& X) const
{
  return this->evalInt(elmInt,fe,X);
}


bool CahnHilliard::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                            const Vec3& X, const Vec3& normal) const
{
  if (!flux)
  {
    std::cerr <<" *** CahnHilliard::evalBou: No flux field."<< std::endl;
    return false;
  }

  double val = normal * (*flux)(X);
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
  if (a) const_cast<CahnHilliard*>(this)->Lnorm = 2;
  return new CahnHilliardNorm(*const_cast<CahnHilliard*>(this),Lnorm,a);
}


CahnHilliard4::CahnHilliard4 (unsigned short int n) : CahnHilliard(n)
{
  scale2nd = 2.0;
}


bool CahnHilliard4::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                             const Vec3& X) const
{
  if (!this->CahnHilliard::evalInt(elmInt,fe,X))
    return false;

  Matrix& A = static_cast<ElmMats&>(elmInt).A.front();
  double s4JxW = pow(smearing,4.0)*fe.detJxW;

  for (size_t i = 1; i <= fe.N.size(); i++)
    for (size_t j = 1; j <= fe.N.size(); j++)
    {
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
  finalOp = Ln == 2 ? ASM::SQRT : ASM::NONE;
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

    if (Lnorm)
    {
      if (C >= 1.0 && (Lnorm == 1 || Lnorm == 2))
        pnorm[k] += fe.detJxW;
      else if (C > 0.0)
      {
        if (Lnorm == 1)
          pnorm[k] += fabs(C)*fe.detJxW; // L1-norm, |c|
        else if (Lnorm == 2)
          pnorm[k] += C*C*fe.detJxW; // L2-norm, |c|
        else if (pnorm[k] == 0.0 || C < pnorm[k])
          pnorm[k] = C; // Smallest-value norm
      }
      k += 2; // Make space for the volume-specific norm |c|/V
    }

    // Dissipated energy, eps_d
    if (C <= 0.0)
      pnorm[k] += 0.25*(Gc/l0)*fe.detJxW;
    else if (C < 1.0)
      pnorm[k] += 0.25*(Gc/l0)*(1.0-C)*(1.0-C)*fe.detJxW;
    pnorm[k++] += Gc*l0*gradC.dot(gradC)*fe.detJxW;

    if (aSol && i == 0)
    {
      // Add final norm group when an analytical solution is provided
      size_t ip = pnorm.size();
      if (aSol->getScalarSecSol() && ip > k+2)
      {
        Vec3 gradCA = (*aSol->getScalarSecSol())(X);
        pnorm[--ip] += (gradCA-gradC).length2()*fe.detJxW;
        pnorm[--ip] += gradCA.length2()*fe.detJxW;
        pnorm[--ip] += gradC*gradC*fe.detJxW;
      }
      if (aSol->getScalarSol() && ip > k+2)
      {
        double CA = (*aSol->getScalarSol())(X);
        pnorm[--ip] += (CA-C)*(CA-C)*fe.detJxW;
        pnorm[--ip] += CA*CA*fe.detJxW;
        pnorm[--ip] += C*C*fe.detJxW;
      }
    }
  }

  return true;
}


bool CahnHilliardNorm::evalIntMx (LocalIntegral& elmInt,
                                  const MxFiniteElement& fe,
                                  const Vec3& X) const
{
  return this->evalInt(elmInt,fe,X);
}


bool CahnHilliardNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (Lnorm == 0) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  size_t nNorm = pnorm.size();
  if (aSol)
  {
    if (aSol->getScalarSol()) nNorm -= 3;
    if (aSol->getScalarSecSol()) nNorm -= 3;
  }

  // Evaluate the volume-specific norm |c|/V
  for (size_t ip = 1; ip < nNorm; ip += 3)
    pnorm[ip+1] = pnorm[ip] / (Lnorm > 0 ? pnorm[0] : 1.0);

  return true;
}


std::string CahnHilliardNorm::getName (size_t i, size_t j,
                                       const char* prefix) const
{
  static const char* errorNorms[] = { "||c^h||_L2", "||c||_L2", "||e^h||_L2",
                                      "||c^h||_H1", "||c||_H1", "||e^h||_H1" };
  std::string name;

  if (aSol && i == this->getNoFields(0))
  {
    if (j <= 6)
      name = errorNorms[j-1];
  }
  else
  {
    if (i == 1 && j == 1)
      return "volume";
    else if (i == 1 && j > 1)
      j --;
    if (Lnorm == 0)
      j += 2;

    if (j == 1)
      name = "|c|";
    else if (j == 2)
      name = "|c|/V";
    else if (j == 3)
      name = "eps_d";
  }

  if (name.empty())
    return this->NormBase::getName(i,j,prefix);
  else if (!prefix)
    return name;

  return std::string(prefix) + " " + name;
}


size_t CahnHilliardNorm::getNoFields (int group) const
{
  size_t nNorm = this->NormBase::getNoFields(0);
  if (group < 1)
    return aSol ? nNorm+1 : nNorm;
  else if (group == (int)(nNorm+1))
  {
    size_t n = 0;
    if (aSol->getScalarSol()) n += 3;
    if (aSol->getScalarSecSol()) n += 3;
    return n;
  }
  else if (group == 1)
    return Lnorm == 0 ? 2 : 4;
  else
    return Lnorm == 0 ? 1 : 3;
}
