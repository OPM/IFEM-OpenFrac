// $Id$
//==============================================================================
//!
//! \file FractureElasticityVoigt.C
//!
//! \date Dec 10 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for elasticity problems with fracture.
//!
//==============================================================================

#include "FractureElasticityVoigt.h"
#include "FiniteElement.h"
#include "MaterialBase.h"
#include "ElmMats.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "Profiler.h"

#ifndef epsZ
//! \brief Zero tolerance for strains.
#define epsZ 1.0e-16
#endif


bool FractureElasticityVoigt::evalStress (double lambda, double mu, double Gc,
                                          const SymmTensor& epsil, double& Phi,
                                          SymmTensor& sigma, Matrix* dSdE) const
{
  PROFILE3("FractureEl::evalStress");

  unsigned short int a = 0, b = 0;

  // Define a Lambda-function to set up the isotropic constitutive matrix
  auto&& setIsotropic = [this,a,b](Matrix& C, double lambda, double mu) mutable
  {
    for (a = 1; a <= C.rows(); a++)
      if (a > nsd)
        C(a,a) = mu;
      else
      {
        C(a,a) = 2.0*mu;
        for (b = 1; b <= nsd; b++)
          C(a,b) += lambda;
      }
  };

  // Define some material constants
  double trEps = epsil.trace();
  double C0 = trEps >= -epsZ ? Gc*lambda : lambda;
  double Cp = Gc*mu;

  if (trEps >= -epsZ && trEps <= epsZ)
  {
    // No strains, stress free configuration
    sigma = Phi = 0.0;
    if (dSdE)
      setIsotropic(*dSdE,C0,Cp);
    return true;
  }

  // Calculate principal strains and the associated directions
  Vec3 eps;
  std::vector<SymmTensor> M(nsd,SymmTensor(nsd));
  {
    PROFILE4("Tensor::principal");
    if (!epsil.principal(eps,M.data()))
      return false;
  }

  // Split the strain tensor into positive and negative parts
  SymmTensor ePos(nsd), eNeg(nsd);
  for (a = 0; a < nsd; a++)
    if (eps[a] > 0.0)
      ePos += eps[a]*M[a];
    else if (eps[a] < 0.0)
      eNeg += eps[a]*M[a];

  // Evaluate the tensile energy
  Phi = mu*(ePos*ePos).trace();
  if (trEps > 0.0) Phi += 0.5*lambda*trEps*trEps;

  // Evaluate the stress tensor
  sigma = C0*trEps;
  sigma += 2.0*mu*(Gc*ePos + eNeg);

#if INT_DEBUG > 4
  std::cout <<"eps_p = "<< eps <<"\n";
  for (a = 0; a < nsd; a++)
    std::cout <<"M("<< 1+a <<") =\n"<< M[a];
  std::cout <<"ePos =\n"<< ePos <<"eNeg =\n"<< eNeg <<"sigma =\n"<< sigma
            <<"Phi = "<< Phi << std::endl;
#endif

  if (!dSdE)
    return true;
  else if (eps[0] == eps[nsd-1])
  {
    // Hydrostatic pressure
    setIsotropic(*dSdE, C0, eps.x > 0.0 ? Cp : mu);
    return true;
  }

  typedef unsigned short int s_ind; // Convenience type definition

  // Define a Lambda-function to calculate (lower triangle of) the matrix Qa
  auto&& getQ = [this](Matrix& Q, const SymmTensor& Ma, double C)
  {
    if (C == 0.0) return;

    auto&& Mult = [Ma](s_ind i, s_ind j, s_ind k, s_ind l)
    {
      return Ma(i,j)*Ma(k,l);
    };

    s_ind i, j, k, l, is, js;

    // Normal stresses and strains
    for (i = 1; i <= nsd; i++)
      for (j = 1; j <= i; j++)
        Q(i,j) += C*Mult(i,i,j,j);

    is = nsd+1;
    for (i = 1; i < nsd; i++)
      for (j = i+1; j <= nsd; j++, is++)
      {
        // Shear stress coupled to normal strain
        for (k = 1; k <= nsd; k++)
          Q(is,k) += C*Mult(i,j,k,k);

        // Shear stress coupled to shear strain
        js = nsd+1;
        for (k = 1; k < nsd; k++)
          for (l = k+1; js <= is; l++, js++)
            Q(is,js) += C*Mult(i,j,k,l);
      }
  };

  // Define a Lambda-function to calculate (lower triangle of) the matrix Gab
  auto&& getG = [this](Matrix& G, const SymmTensor& Ma,
                       const SymmTensor& Mb, double C)
  {
    if (C == 0.0) return;

    auto&& Mult = [Ma,Mb](s_ind i, s_ind j, s_ind k, s_ind l)
    {
      return Ma(i,k)*Mb(j,l) + Ma(i,l)*Mb(j,k) +
             Mb(i,k)*Ma(j,l) + Mb(i,l)*Ma(j,k);
    };

    s_ind i, j, k, l, is, js;

    // Normal stresses and strains
    for (i = 1; i <= nsd; i++)
      for (j = 1; j <= i; j++)
        G(i,j) += C*Mult(i,i,j,j);

    is = nsd+1;
    for (i = 1; i < nsd; i++)
      for (j = i+1; j <= nsd; j++, is++)
      {
        // Shear stress coupled to normal strain
        for (k = 1; k <= nsd; k++)
          G(is,k) += C*Mult(i,j,k,k);

        // Shear stress coupled to shear strain
        js = nsd+1;
        for (k = 1; k < nsd; k++)
          for (l = k+1; js <= is; l++, js++)
            G(is,js) += C*Mult(i,j,k,l);
      }
  };

  // Evaluate the stress tangent assuming Voigt notation and symmetry
  for (a = 1; a <= nsd; a++)
    for (b = 1; b <= a; b++)
      (*dSdE)(a,b) = C0;

  for (a = 0; a < nsd; a++)
  {
    double C1 = eps[a] >= 0.0 ? Cp : mu;
    getQ(*dSdE, M[a], 2.0*C1);
    if (eps[a] != 0.0)
      for (b = 0; b < nsd; b++)
        if (a != b && eps[a] != eps[b])
          getG(*dSdE,M[a],M[b],C1/(1.0-eps[b]/eps[a]));
  }

  // Account for symmetry
  for (b = 2; b <= dSdE->rows(); b++)
    for (a = 1; a < b; a++)
      (*dSdE)(a,b) = (*dSdE)(b,a);

  return true;
}


bool FractureElasticityVoigt::evalInt (LocalIntegral& elmInt,
                                       const FiniteElement& fe,
                                       const Vec3& X) const
{
  PROFILE3("FractureEl::evalInt");

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  size_t nstrc = (nsd+1)*nsd/2;
  Matrix Bmat, dSdE(nstrc,nstrc);
  SymmTensor eps(nsd), sigma(nsd);
  bool lHaveStrains = false;

  if (eKm || eKg || iS || m_mode == SIM::RECOVERY)
  {
    // Evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,0.0,Bmat,eps,eps))
      return false;
    else if (!eps.isZero(1.0e-16))
    {
      lHaveStrains = true;
      // Scale the shear strain components by 0.5 to convert from engineering
      // strains gamma_ij = eps_ij + eps_ji to the tensor components eps_ij
      // which are needed for conistent calculation of principal directions
      for (unsigned short int i = 1; i <= nsd; i++)
        for (unsigned short int j = i+1; j <= nsd; j++)
          eps(i,j) *= 0.5;
    }
#if INT_DEBUG > 3
    std::cout <<"\nFractureElasticity::evalInt(X = "<< X <<")\nBmat ="<< Bmat;
#endif

    // Evaluate the material parameters at this point
    double lambda, mu;
    if (!material->evaluate(lambda,mu,fe,X))
      return false;

    // Evaluate the stress degradation function
    double Gc = this->getStressDegradation(fe.N,elmInt.vec[1]);
#if INT_DEBUG > 3
    std::cout <<"lambda = "<< lambda <<" mu = "<< mu <<" G(c) = "<< Gc <<"\n";
    if (lHaveStrains) std::cout <<"eps =\n"<< eps;
#endif

    // Evaluate the stress state at this point
    if (!this->evalStress(lambda,mu,Gc,eps,myPhi[fe.iGP],sigma,
                          eKm ? &dSdE : nullptr))
      return false;
  }

  if (eKm)
  {
#if INT_DEBUG > 3
    std::cout <<"dSdE ="<< dSdE;
#endif
    // Integrate the material stiffness matrix
    Matrix CB;
    CB.multiply(dSdE,Bmat).multiply(fe.detJxW); // CB = dSdE*B*|J|*w
    elMat.A[eKm-1].multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eKg && lHaveStrains) // Integrate the geometric stiffness matrix
    this->formKG(elMat.A[eKg-1],fe.N,fe.dNdX,0.0,sigma,fe.detJxW);

  if (eM) // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (iS && lHaveStrains)
  {
    // Integrate the internal forces
    sigma *= -fe.detJxW;
    if (!Bmat.multiply(sigma,elMat.b[iS-1],true,true)) // ES -= B^T*sigma
      return false;
  }

  if (eS) // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,fe.detJxW);

  return true;
}


bool FractureElasticityVoigt::evalStress (double lambda, double mu, double Gc,
                                          const SymmTensor& epsilon,
                                          double& Phi, SymmTensor& sigma) const
{
  return this->evalStress(lambda,mu,Gc,epsilon,Phi,sigma,nullptr);
}
