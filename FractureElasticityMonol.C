// $Id$
//==============================================================================
//!
//! \file FractureElasticityMonol.C
//!
//! \date Jul 10 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementation for monolithic fracture elasticity problems.
//!
//==============================================================================

#include "FractureElasticityMonol.h"
#include "FiniteElement.h"
#include "BlockElmMats.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


FractureElasticityMonol::FractureElasticityMonol (unsigned short int n, int ord)
  : FractureElasticityVoigt(n)
{
  Gc = smearing = gammaInv = 1.0;
  crtol = 0.0;
  use4th = ord == 4;
  npv = nsd + 1;
}


bool FractureElasticityMonol::parse (const TiXmlElement* elem)
{
  const char* value = utl::getValue(elem,"Gc");
  if (value)
    Gc = atof(value);
  else if ((value = utl::getValue(elem,"smearing")))
    smearing = atof(value);
  else if ((value = utl::getValue(elem,"penalty_factor"))) {
    gammaInv = 1.0/atof(value);
    utl::getAttribute(elem,"threshold",crtol);
  }
  else
    return this->FractureElasticityVoigt::parse(elem);

  return true;
}


void FractureElasticityMonol::printLog () const
{
  this->FractureElasticityVoigt::printLog();

  IFEM::cout <<"\tCritical fracture energy density: "<< Gc
             <<"\n\tSmearing factor: "<< smearing
             <<"\n\tUsing monolithic coupling to "
             << (use4th ? "fourth" : "second")
             <<"-order phase field."<< std::endl;
  IFEM::cout <<"\tPenalty parameter: "<< 1.0/gammaInv
             <<" (threshold value: "<< crtol <<")"<< std::endl;
}


void FractureElasticityMonol::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;

  eM = eKm = eKg = eKc = eAcc = 0;
  eS = iS = eBc = 0;

  switch (mode)
    {
    case SIM::STATIC:
      eKm  = eKg = 2;
      eS   = iS  = 2;
      eAcc = 3;
      eBc  = 3;
      eKc  = 4;
      if (intPrm[3] > 0.0)
        eKg = 0; // Linear analysis, no geometric stiffness
      primsol.resize(nSV);
      break;

    case SIM::RHS_ONLY:
      eS  = iS = 2;
      eBc = 3;

    case SIM::RECOVERY:
      primsol.resize(1);
      break;

    default:
      primsol.clear();
    }
}


LocalIntegral* FractureElasticityMonol::getLocalIntegral (size_t nen, size_t,
                                                          bool neumann) const
{
  BlockElmMats* result = new BlockElmMats(2);

  switch (m_mode)
  {
    case SIM::STATIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann ? 0 : 4, 3);
      break;

    case SIM::RHS_ONLY:
      result->resize(neumann ? 0 : 4, 3);

    case SIM::RECOVERY:
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      ; // Other solution modes not supported
  }

  result->redim(1,nen,nsd);
  result->redim(2,nen,1);
  result->redimOffDiag(3);
  result->redimNewtonMat();

  return result;
}


bool FractureElasticityMonol::initElement (const std::vector<int>& MNPC,
                                           LocalIntegral& elmInt)
{
  if (mySol.empty())
  {
    std::cerr <<" *** FractureElasticityMonol::initElement:"
              <<" No primary solution vectors."<< std::endl;
    return false;
  }
  else if (m_mode == SIM::DYNAMIC)
  {
    std::cerr <<" *** FractureElasticityMonol::initElement:"
              <<" Monolithic coupling is not available"
              <<" for dynamic simulation."<< std::endl;
    return false;
  }

  size_t nsol = 1 + eC;
  if (elmInt.vec.size() < nsol)
    elmInt.vec.resize(nsol);

  return this->getSolution(elmInt.vec,MNPC);
}


bool FractureElasticityMonol::evalInt (LocalIntegral& elmInt,
                                       const FiniteElement& fe,
                                       const Vec3& X) const
{
  if (!this->FractureElasticityVoigt::evalInt(elmInt,fe,X))
    return false;
  else if (!eAcc && !eBc)
    return true;

  double C = fe.N.dot(elmInt.vec[eC]);
  double ddGc = this->getStressDegradation(fe.N,elmInt.vec,2);
  double scale = ddGc*myPhi[fe.iGP] + 0.5*Gc/smearing;
  if (C < crtol) scale += gammaInv;
  double s1JxW = scale*fe.detJxW;
  double s2JxW = (use4th ? 1.0 : 2.0)*Gc*smearing*fe.detJxW;

  if (eAcc)
  {
    Matrix& A = static_cast<ElmMats&>(elmInt).A[eAcc-1];

    for (size_t i = 1; i <= fe.N.size(); i++)
      for (size_t j = 1; j <= fe.N.size(); j++)
      {
        double grad = 0.0;
        for (size_t k = 1; k <= nsd; k++)
          grad += fe.dNdX(i,k)*fe.dNdX(j,k);
        A(i,j) += fe.N(i)*fe.N(j)*s1JxW + grad*s2JxW;
      }

    if (use4th)
    {
      // Forth-order phase field
      double s4JxW = 0.5*Gc*pow(smearing,3.0)*fe.detJxW;

      for (size_t i = 1; i <= fe.N.size(); i++)
        for (size_t j = 1; j <= fe.N.size(); j++)
        {
          double grad = 0.0;
          for (unsigned short int k = 1; k <= nsd; k++)
            grad += fe.d2NdX2(i,k,k)*fe.d2NdX2(j,k,k);
          A(i,j) += grad*s4JxW;
        }
    }
  }

  if (eBc)
  {
    Vector& R = static_cast<ElmMats&>(elmInt).b[eBc-1];

    double dGc = this->getStressDegradation(fe.N,elmInt.vec,1);
    scale = dGc*myPhi[fe.iGP] - 0.5*(1.0-C)*Gc/smearing;
    if (C < crtol) scale += C*gammaInv;
    R.add(fe.N,-scale*fe.detJxW); // R -= N*scale*detJxW

    Vector gradC; // Compute the phase field gradient gradC = dNdX^t*eC
    if (!fe.dNdX.multiply(elmInt.vec[eC],gradC,true))
      return false;

    gradC *= -s2JxW;
    return fe.dNdX.multiply(gradC,R,false,true); // R -= dNdX*gradC*s2JxW
  }

  return true;
}


bool FractureElasticityMonol::evalSol (Vector& s,
                                       const FiniteElement& fe, const Vec3& X,
                                       const std::vector<int>& MNPC) const
{
  Vectors eV(1+eC);
  if (!mySol.empty() && this->getSolution(eV,MNPC))
    return this->evalSol2(s,eV,fe,X);
  else
    return false;
}


bool FractureElasticityMonol::getSolution (Vectors& eV,
                                           const std::vector<int>& MNPC) const
{
  // Extract element displacements and phase field from the monolithic solution
  Matrix temp(npv,MNPC.size());
  int ierr = utl::gather(MNPC,npv,mySol.front(),temp);
  if (ierr > 0)
  {
    std::cerr <<" *** FractureElasticityMonol::getSolution: Detected "
              << ierr <<" node numbers out of range."<< std::endl;
    return false;
  }

  eV[eC]     = temp.getRow(npv);
  eV.front() = temp.expandRows(-1);
#if INT_DEBUG > 2
  std::cout <<"Element displacement vector:"<< eV.front()
            <<"Element phase field vector:"<< eV[eC];
#endif
  return true;
}


size_t FractureElasticityMonol::getNoFields (int fld) const
{
  return fld < 2 ? npv : this->FractureElasticityVoigt::getNoFields(fld);
}


std::string FractureElasticityMonol::getField1Name (size_t i,
                                                    const char* prefix) const
{
  if (i < nsd)
    return this->FractureElasticityVoigt::getField1Name(i,prefix);
  else if (i > nsd)
    return this->FractureElasticityVoigt::getField1Name(4,prefix);

  return prefix ? prefix + std::string(" phase") : std::string("phase");
}
