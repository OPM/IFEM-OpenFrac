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
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "Profiler.h"

#ifndef epsZ
//! \brief Zero tolerance for strains.
#define epsZ 1.0e-16
#endif


LocalIntegral* FractureElasticityVoigt::getLocalIntegral (size_t nen, size_t,
                                                          bool neumann) const
{
  LocalIntegral* li = this->FractureElasticity::getLocalIntegral(nen,0,neumann);
  if (m_mode >= SIM::RHS_ONLY && !neumann)
    static_cast<ElmMats*>(li)->c.resize(1); // Total strain energy
  return li;
}


bool FractureElasticityVoigt::evalStress (double lambda, double mu, double Gc,
                                          const SymmTensor& epsil, double* Phi,
                                          SymmTensor* sigma, Matrix* dSdE,
                                          bool postProc, bool printElm) const
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
    Phi[0] = 0.0;
    if (postProc)
      Phi[1] = Phi[2] = Phi[3] = 0.0;
    if (sigma)
      *sigma = 0.0;
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

  if (sigma)
  {
    // Evaluate the stress tensor
    *sigma = C0*trEps;
    *sigma += 2.0*mu*(Gc*ePos + eNeg);
  }

  // Evaluate the tensile energy
  Phi[0] = mu*(ePos*ePos).trace();
  if (trEps > 0.0) Phi[0] += 0.5*lambda*trEps*trEps;
  // Evaluate the compressive energy
  Phi[1] = mu*(eNeg*eNeg).trace();
  if (trEps < 0.0) Phi[1] += 0.5*lambda*trEps*trEps;
  // Evaluate the total strain energy
  Phi[2] = Gc*Phi[0] + Phi[1];
  if (postProc) // Evaluate the bulk energy
    Phi[3] = Gc*(Phi[0] + Phi[1]);
  else if (sigmaC > 0.0) // Evaluate the crack driving function
    Phi[0] = this->MieheCrit56(eps,lambda,mu);

#if INT_DEBUG > 4
  std::cout <<"eps_p = "<< eps <<"\n";
  for (a = 0; a < nsd; a++)
    std::cout <<"M("<< 1+a <<") =\n"<< M[a];
  std::cout <<"ePos =\n"<< ePos <<"eNeg =\n"<< eNeg;
  if (sigma) std::cout <<"sigma =\n"<< *sigma;
  std::cout <<"Phi = "<< Phi[0] <<" "<< Phi[1] <<" "<< Phi[2];
  if (postProc) std::cout <<" "<< Phi[3];
  std::cout << std::endl;
#else
  if (printElm)
  {
    std::cout <<"g(c) = "<< Gc
              <<"\nepsilon =\n"<< epsil <<"eps_p = "<< eps
              <<"\nePos =\n"<< ePos <<"eNeg =\n"<< eNeg;
    if (sigma) std::cout <<"sigma =\n"<< *sigma;
    std::cout <<"Phi = "<< Phi[0] <<" "<< Phi[1] <<" "<< Phi[2];
    if (postProc) std::cout <<" "<< Phi[3];
    std::cout << std::endl;
  }
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
  double U = 0.0;

  if (eKm || eKg || iS || m_mode == SIM::RECOVERY)
  {
    // Evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,0.0,Bmat,eps,eps))
      return false;
    else if (!eps.isZero(1.0e-16))
      lHaveStrains = true;

#if INT_DEBUG > 3
    std::cout <<"\nFractureElasticity::evalInt(X = "<< X <<")\nBmat ="<< Bmat;
#endif

    // Evaluate the stress degradation function
    double Gc = this->getStressDegradation(fe.N,elmInt.vec);
    if (tSplit < 0.0 || static_cast<const Vec4&>(X).t < tSplit)
    {
      // Evaluate the constitutive matrix and the stress tensor at this point
      if (!material->evaluate(dSdE,sigma,U,fe,X,eps,eps,3))
        return false;

      // Degrade the stresses and strain energy isotropically
      dSdE *= Gc;
      sigma *= Gc;
      myPhi[fe.iGP] = (U *= Gc);

#if INT_DEBUG > 3
      std::cout <<"G(c) = "<< Gc <<"\n";
      if (lHaveStrains)
        std::cout <<"eps =\n"<< eps <<"sigma =\n"<< sigma
                  <<"Phi = "<< myPhi[fe.iGP] << std::endl;
#endif
    }
    else
    {
      // Evaluate the material parameters at this point
      double lambda, mu;
      if (!material->evaluate(lambda,mu,fe,X))
        return false;

      // Scale the shear strain components by 0.5 to convert from engineering
      // strains gamma_ij = eps_ij + eps_ji to the tensor components eps_ij
      // which are needed for consistent calculation of principal directions
      for (unsigned short int i = 1; i <= nsd; i++)
        for (unsigned short int j = i+1; j <= nsd; j++)
          eps(i,j) *= 0.5;

#if INT_DEBUG > 3
      std::cout <<"lambda = "<< lambda <<" mu = "<< mu <<" G(c) = "<< Gc <<"\n";
      if (lHaveStrains) std::cout <<"eps =\n"<< eps;
#endif

      // Evaluate the stress state at this point, with degraded tensile part
      double Phi[3];
      if (!this->evalStress(lambda,mu,Gc,eps,Phi,&sigma,
                            eKm ? &dSdE : nullptr))
        return false;

      myPhi[fe.iGP] = Phi[0];
      U = Phi[2];
    }
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

  if (eS)
  {
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,fe.detJxW);
    // Integrate the load vector due to internal crack pressure
    this->formCrackForce(elMat.b[eS-1],elMat.vec,fe,X);
  }

  if (lHaveStrains && !elMat.c.empty())
    // Integrate the total strain energy
    elMat.c.front() += U*fe.detJxW;

  return true;
}


bool FractureElasticityVoigt::evalStress (double lambda, double mu, double Gc,
                                          const SymmTensor& epsilon,
                                          double* Phi, SymmTensor& sigma) const
{
  return this->evalStress(lambda,mu,Gc,epsilon,Phi,&sigma,nullptr,true);
}


NormBase* FractureElasticityVoigt::getNormIntegrand (AnaSol*) const
{
  return new FractureElasticNorm(*const_cast<FractureElasticityVoigt*>(this));
}


int FractureElasticNorm::dbgElm = 0;


FractureElasticNorm::FractureElasticNorm (FractureElasticityVoigt& p)
  : ElasticityNorm(p)
{
  finalOp = ASM::NONE;
}


size_t FractureElasticNorm::getNoFields (int group) const
{
  return group < 1 ? 1 : 5;
}


std::string FractureElasticNorm::getName (size_t, size_t j, const char*) const
{
  if (j == 1)
    return "Elastic strain energy";
  else if (j == 2)
    return "External energy";
  else if (j == 3)
    return "Tensile energy";
  else if (j == 4)
    return "Compressive energy";
  else if (j == 5)
    return "Bulk energy";
  else
    return this->NormBase::getName(1,j);
}


bool FractureElasticNorm::evalInt (LocalIntegral& elmInt,
                                   const FiniteElement& fe,
                                   const Vec3& X) const
{
  FractureElasticityVoigt& p = static_cast<FractureElasticityVoigt&>(myProblem);
  ElmNorm&             pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the symmetric strain tensor, eps
  Matrix Bmat;
  SymmTensor eps(p.getNoSpaceDim());
  if (!p.kinematics(elmInt.vec.front(),fe.N,fe.dNdX,0.0,Bmat,eps,eps))
    return false;
  else if (!eps.isZero(1.0e-16))
    // Scale the shear strain components by 0.5 to convert from engineering
    // strains gamma_ij = eps_ij + eps_ji to the tensor components eps_ij
    // which are needed for consistent calculation of principal directions
    for (unsigned short int i = 1; i <= eps.dim(); i++)
      for (unsigned short int j = i+1; j <= eps.dim(); j++)
        eps(i,j) *= 0.5;

  bool printElm = fe.iel == dbgElm;
  if (printElm)
    std::cout <<"\nFractureElasticNorm::evalInt: iel,ip,X = "
              << fe.iel <<" "<< fe.iGP <<" "<< X << std::endl;

  // Evaluate the strain energy at this point
  double Phi[4];
  double Gc = p.getStressDegradation(fe.N,elmInt.vec);
  if (p.tSplit < 0.0 || static_cast<const Vec4&>(X).t < p.tSplit)
  {
    // Evaluate the strain energy density at this point
    SymmTensor sigma(eps.dim());
    if (!p.material->evaluate(Bmat,sigma,Phi[2],fe,X,eps,eps,3))
      return false;
    Phi[2] *= Gc; // Isotropic scaling
    Phi[0] = Phi[1] = Phi[3] = 0.0;
  }
  else
  {
    // Evaluate the material parameters at this point
    double lambda, mu;
    if (!p.material->evaluate(lambda,mu,fe,X))
      return false;
    // Evaluate the tensile-degraded strain energy
    if (!p.evalStress(lambda,mu,Gc,eps,Phi,nullptr,nullptr,true,printElm))
      return false;
  }

  // Integrate the total elastic energy
  pnorm[0] += Phi[2]*fe.detJxW;

  if (p.haveLoads())
  {
    // Evaluate the body load
    Vec3 f = p.getBodyforce(X);
    // Evaluate the displacement field
    Vec3 u = p.evalSol(pnorm.vec.front(),fe.N);
    // Integrate the external energy (f,u^h)
    pnorm[1] += f*u*fe.detJxW;
  }

  // Integrate the tensile and compressive energies
  pnorm[2] += Phi[0]*fe.detJxW;
  pnorm[3] += Phi[1]*fe.detJxW;
  // Integrate the bulk energy
  pnorm[4] += Phi[3]*fe.detJxW;

  return true;
}
