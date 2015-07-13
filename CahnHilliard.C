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
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "VTF.h"
#include "IFEM.h"


CahnHilliard::CahnHilliard(unsigned short int n) :
  nsd(n), smearFactor(1.0), maxCrack(1e-3), stabk(0.0),
  initial_crack(nullptr), mat(nullptr), tensile(nullptr)
{
  npv = 1; // One primary unknown per node (scalar equation)

  primsol.resize(1);
}


LocalIntegral* CahnHilliard::getLocalIntegral (size_t nen, size_t,
                                          bool neumann) const
{
  ElmMats* result = new ElmMats();
  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->resize(neumann ? 0 : 1, 1);

  result->redim(nen);
  return result;
}


bool CahnHilliard::evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                           const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  double Gc = mat->getFractureEnergy(X);

  if (initial_crack) {
    double dist = (*initial_crack)(X);
    if (dist < smearFactor)
      historyField(fe.iGP+1) = Gc/(4.0*smearFactor)*(1.0/maxCrack-1.0)*(1.0-dist/smearFactor);
  }

  // update history field
  historyField(fe.iGP+1) = std::max(historyField(fe.iGP+1), tensile?(*tensile)(fe.iGP+1):0.0);

  double scale = 4.0*smearFactor*(1.0-stabk)*historyField(fe.iGP+1)/Gc;

  for (size_t i=1;i<=fe.N.size();++i) {
    for (size_t j=1;j<=fe.N.size();++j) {
      elMat.A[0](i,j) += (scale+1.0)*fe.N(i)*fe.N(j)*fe.detJxW;
      double grad=0.0;
      for (size_t k=1;k<=nsd;++k)
        grad += fe.dNdX(i,k)*fe.dNdX(j,k);
      elMat.A[0](i,j) += 4.0*smearFactor*smearFactor*grad*fe.detJxW;
    }
  }

  elMat.b[0].add(fe.N,fe.detJxW);

  return true;
}


void CahnHilliard::initIntegration (size_t nIp, size_t nBp)
{
  historyField.resize(nIp);
}


std::string CahnHilliard::getField1Name(size_t, const char* prefix) const
{
  if (!prefix) return "c";

  return std::string(prefix) + " c";
}


void CahnHilliard::printLog()
{
  IFEM::cout << "Smearing factor: " << smearFactor << std::endl;
  IFEM::cout << "Max value in crack: " << maxCrack << std::endl;
}


bool CahnHilliard::evalSol(Vector& s, const FiniteElement& fe,
                           const Vec3& X, const std::vector<int>& MNPC) const
{
  s.resize(1);
  Vector tmp;
  if (utl::gather(MNPC,1,primsol.front(),tmp))
    return false;

  s(1) = fe.N.dot(tmp);
  if (s(1) < maxCrack)
    s(1) = 0.0;
  if (s(1) > 1.0)
    s(1) = 1.0;

  return true;
}
