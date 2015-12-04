// $Id$
//==============================================================================
//!
//! \file FractureElasticity.C
//!
//! \date Oct 27 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for elasticity problems with fracture.
//!
//==============================================================================

#include "FractureElasticity.h"
#include "FiniteElement.h"
#include "MaterialBase.h"
#include "ElmMats.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Tensor4.h"
#include "Tensor.h"

#ifndef epsZ
//! \brief Zero tolerance for strains.
#define epsZ 1.0e-16
#endif


FractureElasticity::FractureElasticity (unsigned short int n) : Elasticity(n)
{
  alpha = 0.0;
  myPhi = nullptr;
  this->registerVector("phasefield",&myCVec);
}


void FractureElasticity::initIntegration (size_t nGp, size_t)
{
  // Initialize internal tensile energy buffer
  delete[] myPhi;
  myPhi = new double[nGp];
  memset(myPhi,0,nGp*sizeof(double));
}


bool FractureElasticity::initElement (const std::vector<int>& MNPC,
                                      LocalIntegral& elmInt)
{
  if (primsol.empty())
  {
    std::cerr <<" *** FractureElasticity::initElement:"
              <<" No primary solution vectors."<< std::endl;
    return false;
  }

  int ierr = 0;
  elmInt.vec.resize(primsol.size()+1);

  // Extract displacement vector for this element
  if (!primsol.front().empty())
    ierr = utl::gather(MNPC,npv,primsol.front(),elmInt.vec.front());

  // Extract crack phase field vector for this element
  if (ierr == 0 && !myCVec.empty())
    ierr = utl::gather(MNPC,1,myCVec,elmInt.vec[1]);

  // Extract velocity and acceleration vectors for this element
  for (size_t i = 1; i < primsol.size() && ierr == 0; i++)
    if (!primsol[i].empty())
      ierr = utl::gather(MNPC,npv,primsol[i],elmInt.vec[1+i]);

  if (ierr == 0) return true;

  std::cerr <<" *** FractureElasticity::initElement: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


void FractureElasticity::getG (const Tensor& Ma, const Tensor& Mb,
                               Tensor4& Gab, double eps) const
{
  unsigned short int i, j, k, l;
  for (i = 1; i <= nsd; i++)
    for (j = 1; j <= nsd; j++)
      for (k = 1; k <= nsd; k++)
        for (l = 1; l <= nsd; l++)
          Gab(i,j,k,l) += eps*(Ma(i,k)*Mb(j,l) + Ma(i,l)*Mb(j,k) +
                               Mb(i,k)*Ma(j,l) + Mb(i,l)*Ma(j,k));
}


void FractureElasticity::getQ (const Tensor& Ma, Tensor4& Qa, double C) const
{
  unsigned short int i, j, k, l;
  for (i = 1; i <= nsd; i++)
    for (j = 1; j <= nsd; j++)
      for (k = 1; k <= nsd; k++)
        for (l = 1; l <= nsd; l++)
          Qa(i,j,k,l) += C*Ma(i,j)*Ma(k,l);
}


bool FractureElasticity::evalStress (double lambda, double mu, double Gc,
                                     const SymmTensor& epsilon, double& Phi,
                                     SymmTensor& sigma, Tensor4& dSdE) const
{
  // Define some material constants
  double trEps = epsilon.trace();
  double C0 = trEps >= -epsZ ? Gc*lambda : lambda;
  double Cm = 2.0*mu;
  double Cp = Gc*Cm;
  unsigned short int a, b;

  // Set up the stress tangent (4th order tensor)
  dSdE = Tensor4(nsd,C0,true);
  if (trEps >= -epsZ && trEps <= epsZ)
  {
    // No strains, stress free configuration
    sigma = Phi = 0.0;
    for (a = 1; a <= nsd; a++)
      for (b = 1; b <= nsd; b++)
        dSdE(a,b,a,b) += Cp;
    return true;
  }

  // Calculate principal strains and the associated directions
  Vec3 eps;
  std::vector<SymmTensor> M(nsd,SymmTensor(nsd));
  if (!epsilon.principal(eps,M.data()))
    return false;

  // Split the strain tensor into positive and negative part
  SymmTensor ePos(nsd), eNeg(nsd);
  for (a = 0; a < nsd; a++)
    if (eps[a] > 0.0)
      ePos += eps[a]*M[a];
    else if (eps[a] < 0.0)
      eNeg += eps[a]*M[a];

  // Evaluate the tensile energy
  Phi = 0.5*(C0*trEps*trEps + Cp*(ePos*ePos).trace() + Cm*(eNeg*eNeg).trace());

  // Evaluate the stress tensor
  sigma = C0*trEps;
  sigma += Cp*ePos + Cm*eNeg;

#if INT_DEBUG > 4
  std::cout <<"eps_p = "<< eps <<"\n";
  for (a = 0; a < nsd; a++)
    std::cout <<"M("<< 1+a <<") =\n"<< M[a];
  std::cout <<"ePos =\n"<< ePos <<"eNeg =\n"<< eNeg <<"sigma =\n"<< sigma
            <<"Phi = "<< Phi << std::endl;
#endif

  if (eps[0] == eps[nsd-1])
  {
    // Hydrostatic pressure
    for (a = 1; a <= nsd; a++)
      for (b = 1; b <= nsd; b++)
        dSdE(a,b,a,b) += eps.x > 0.0 ? Cp : Cm;
    return true;
  }

  // Evaluate the stress tangent (4th order tensor)
  for (a = 0; a < nsd; a++)
  {
    if (eps[a] >= 0.0)
      this->getQ(M[a],dSdE,Cp);
    if (eps[a] <= 0.0)
      this->getQ(M[a],dSdE,Cm);
    for (b = 0; b < nsd; b++)
      if (a != b && eps[a] != eps[b])
      {
        double epsRel = 0.5*eps[a] / (eps[a]-eps[b]);
        if (eps[a] >= 0.0)
          this->getG(M[a],M[b],dSdE,Cp*epsRel);
        if (eps[a] <= 0.0)
          this->getG(M[a],M[b],dSdE,Cm*epsRel);
      }
  }

  return true;
}


bool FractureElasticity::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  Matrix Bmat;
  SymmTensor eps(nsd), sigma(nsd);
  bool lHaveStrains = false;

  if (eKm || eKg)
  {
    // Evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,0.0,Bmat,eps,eps))
      return false;
    else if (!eps.isZero(1.0e-16))
      lHaveStrains = true;
#if INT_DEBUG > 3
    std::cout <<"\nFractureElasticity::evalInt(X = "<< X <<")\nBmat ="<< Bmat;
#endif

    // Evaluate the material parameters at this point
    double lambda, mu;
    if (!material->evaluate(lambda,mu,fe,X))
      return false;

    // Evaluate the crack phase field function
    double cc = elmInt.vec[1].empty() ? 1.0 : fe.N.dot(elmInt.vec[1]);
    double Gc = (1.0-alpha)*cc*cc + alpha;
#if INT_DEBUG > 3
    std::cout <<"lambda = "<< lambda <<" mu = "<< mu <<" G(c) = "<< Gc <<"\n";
    if (lHaveStrains) std::cout <<"eps =\n"<< eps;
#endif

    // Evaluate the stress state at this point
    Tensor4 dSdE(nsd);
    if (!this->evalStress(lambda,mu,Gc,eps,myPhi[fe.iGP],sigma,dSdE))
      return false;
#if INT_DEBUG > 3
    std::cout <<"dSdE ="<< dSdE;
#endif

    // Integrate the material stiffness matrix
    Matrix& Km = elMat.A[eKm-1];
    size_t a, b, nstrc = Bmat.rows();
    unsigned short int i, j, k, l;
    for (b = 1; b <= Km.cols(); b++)
    {
      const double* Bb = Bmat.ptr(b-1);
      SymmTensor dEpsBj(RealArray(Bb,Bb+nstrc));
      for (i = 1; i <= nsd; i++)
        for (j = i+1; j <= nsd; j++)
          dEpsBj(i,j) *= 0.5;

      Matrix Ctmp(nsd,nsd);
      for (i = 1; i <= nsd; i++)
        for (j = 1; j <= nsd; j++)
          for (k = 1; k <= nsd; k++)
            for (l = 1; l <= nsd; l++)
              Ctmp(i,j) += dSdE(i,j,k,l)*dEpsBj(k,l);

      for (a = 1; a <= Km.rows(); a++)
      {
        const double* Ba = Bmat.ptr(a-1);
        SymmTensor dEpsAi(RealArray(Ba,Ba+nstrc));
        for (i = 1; i <= nsd; i++)
          for (j = i+1; j <= nsd; j++)
            dEpsAi(i,j) *= 0.5;
        for (i = 1; i <= nsd; i++)
          for (j = 1; j <= nsd; j++)
            Km(a,b) += dEpsAi(i,j)*Ctmp(i,j)*fe.detJxW;
      }
    }
  }

  if (eKg && lHaveStrains) // Integrate the geometric stiffness matrix
    this->formKG(elMat.A[eKg-1],fe.N,fe.dNdX,0.0,sigma,fe.detJxW);

  if (eM) // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (eS) // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,fe.detJxW);

  return true;
}


bool FractureElasticity::evalBou (LocalIntegral& elmInt,
                                  const FiniteElement& fe, const Vec3& X,
                                  const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr <<" *** FractureElasticity::evalBou: No tractions."<< std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** FractureElasticity::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec3 T = this->getTraction(X,normal);

  // Integrate the force vector
  Vector& ES = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*fe.N(a)*fe.detJxW;

  return true;
}
