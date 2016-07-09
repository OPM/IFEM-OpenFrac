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
  Gc = smearing = 1.0;
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

  result->redim(eKm-1,nen,nsd);
  result->redim(eAcc-1,nen,1);
  result->redimOffDiag(eKc-1);
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

  // Extract the element level solution vectors from the monolithic solution
  Matrix temp(npv,MNPC.size());
  int ierr = utl::gather(MNPC,npv,primsol.front(),temp);
  elmInt.vec[eC] = temp.getRow(npv);
  elmInt.vec.front() = temp.expandRows(-1);

  if (ierr == 0) return true;

  std::cerr <<" *** FractureElasticityMonol::initElement: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


bool FractureElasticityMonol::evalInt (LocalIntegral& elmInt,
                                       const FiniteElement& fe,
                                       const Vec3& X) const
{
  if (!this->FractureElasticityVoigt::evalInt(elmInt,fe,X))
    return false;

  if (eAcc)
  {
    Matrix& A = static_cast<ElmMats&>(elmInt).A[eAcc-1];

    double PhiPl = myPhi[fe.iGP];
    double scale = 1.0 + 4.0*smearing*PhiPl/Gc;
    double s1JxW = scale*fe.detJxW;
    double s2JxW = (use4th ? 2.0 : 4.0)*smearing*smearing*fe.detJxW;

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
      double s4JxW = pow(smearing,4.0)*fe.detJxW;

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
    static_cast<ElmMats&>(elmInt).b[eBc-1].add(fe.N,fe.detJxW);

  return true;
}


bool FractureElasticityMonol::evalSol (Vector& s,
                                       const FiniteElement& fe, const Vec3& X,
                                       const std::vector<int>& MNPC) const
{
  // Extract element displacements and phase field from the monolithic solution
  Vectors eV(1+eC);
  Matrix temp(npv,MNPC.size());
  int ierr = utl::gather(MNPC,npv,mySol.front(),temp);
  eV[eC] = temp.getRow(npv);
  eV.front() = temp.expandRows(-1);

  if (ierr > 0)
  {
    std::cerr <<" *** FractureElasticityMonol::evalSol: Detected "<< ierr
              <<" node numbers out of range."<< std::endl;
    return false;
  }

  return this->evalSol2(s,eV,fe,X);
}
