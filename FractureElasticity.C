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
#include "Profiler.h"

#ifndef epsZ
//! \brief Zero tolerance for strains.
#define epsZ 1.0e-16
#endif


FractureElasticity::FractureElasticity (unsigned short int n) : Elasticity(n)
{
  alpha = 0.0;
  this->registerVector("phasefield",&myCVec);
}


void FractureElasticity::initIntegration (size_t nGp, size_t)
{
  // Initialize internal tensile energy buffer
  myPhi.resize(nGp);
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


bool FractureElasticity::evalStress (double lambda, double mu, double Gc,
                                     const SymmTensor& epsilon, double& Phi,
                                     SymmTensor& sigma, Tensor4* dSdE) const
{
  PROFILE3("FractureEl::evalStress");

  unsigned short int a = 0, b = 0;

  // Define a Lambda-function to set up the isotropic constitutive tensor
  auto&& setIsotropic = [this,a,b](Tensor4& C, double mu) mutable
  {
    for (a = 1; a <= nsd; a++)
      for (b = 1; b <= nsd; b++)
      {
        C(a,b,a,b) += mu;
        C(a,b,b,a) += mu;
      }
  };

  // Define some material constants
  double trEps = epsilon.trace();
  double C0 = trEps >= -epsZ ? Gc*lambda : lambda;
  double Cp = Gc*mu;

  // Set up the stress tangent (4th order tensor)
  if (dSdE)
    *dSdE = Tensor4(nsd,C0,true);
  if (trEps >= -epsZ && trEps <= epsZ)
  {
    // No strains, stress free configuration
    sigma = Phi = 0.0;
    if (dSdE)
      setIsotropic(*dSdE,Cp);
    return true;
  }

  // Calculate principal strains and the associated directions
  Vec3 eps;
  std::vector<SymmTensor> M(nsd,SymmTensor(nsd));
  {
    PROFILE4("Tensor::principal");
    if (!epsilon.principal(eps,M.data()))
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
  Phi = 0.5*C0*trEps*trEps + mu*(Gc*(ePos*ePos).trace() + (eNeg*eNeg).trace());

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

  if (!dSdE) return true;

  if (eps[0] == eps[nsd-1])
  {
    // Hydrostatic pressure
    setIsotropic(*dSdE, eps.x > 0.0 ? Cp : mu);
    return true;
  }

  // Define a Lambda-function to calculate the tensor Qa
  auto&& getQ = [this](Tensor4& Qa, const Tensor& Ma, double C)
  {
    if (C == 0.0) return;

    unsigned short int i, j, k, l;
    for (i = 1; i <= nsd; i++)
      for (j = 1; j <= nsd; j++)
        for (k = 1; k <= nsd; k++)
          for (l = 1; l <= nsd; l++)
            Qa(i,j,k,l) += C*Ma(i,j)*Ma(k,l);
  };

  // Define a Lambda-function to calculate the tensor Gab
  auto&& getG = [this](Tensor4& Gab, const Tensor& Ma,
                       const Tensor& Mb, double eps)
  {
    if (eps == 0.0) return;

    unsigned short int i, j, k, l;
    for (i = 1; i <= nsd; i++)
      for (j = 1; j <= nsd; j++)
        for (k = 1; k <= nsd; k++)
          for (l = 1; l <= nsd; l++)
            Gab(i,j,k,l) += eps*(Ma(i,k)*Mb(j,l) + Ma(i,l)*Mb(j,k) +
                                 Mb(i,k)*Ma(j,l) + Mb(i,l)*Ma(j,k));
  };

  // Evaluate the stress tangent (4th order tensor)
  for (a = 0; a < nsd; a++)
  {
    double C1 = eps[a] >= 0.0 ? Cp : mu;
    getQ(*dSdE, M[a], 2.0*C1);
    if (eps[a] != 0.0)
      for (b = 0; b < nsd; b++)
        if (a != b && eps[a] != eps[b])
          getG(*dSdE, M[a], M[b], C1/(1.0-eps[b]/eps[a]));
  }

  return true;
}


double FractureElasticity::getStressDegradation (const Vector& N,
                                                 const Vector& eC) const
{
  // Evaluate the crack phase field function, c(X)
  double c = eC.empty() ? 1.0 : N.dot(eC);
  // Evaluate the stress degradation function, g(c), ignoring negative values
  return c > 0.0 ? (1.0-alpha)*c*c + alpha : alpha;
}


bool FractureElasticity::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe, const Vec3& X) const
{
  PROFILE3("FractureEl::evalInt");

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  Matrix Bmat;
  Tensor4 dSdE(nsd);
  SymmTensor eps(nsd), sigma(nsd);
  bool lHaveStrains = false;

  if (eKm || eKg || iS)
  {
    // Evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,0.0,Bmat,eps,eps))
      return false;
    else if (!eps.isZero(1.0e-16))
    {
      lHaveStrains = true;
      for (unsigned short int i = 1; i <= nsd; i++)
        for (unsigned short int j = i+1; j <= nsd; j++)
          eps(i,j) *= 0.5; // Using tensor formulation
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

  // Define a Lambda-function to extract a symmetric tensor from the B-matrix
  auto&& getDepsDu = [this,&Bmat](SymmTensor& dEpsDu, size_t a)
  {
    const double* Ba = Bmat.ptr(a-1);
    dEpsDu = RealArray(Ba,Ba+Bmat.rows());
    for (unsigned short int i = 1; i <= nsd; i++)
      for (unsigned short int j = i+1; j <= nsd; j++)
        dEpsDu(i,j) *= 0.5;
  };

  if (eKm)
  {
#if INT_DEBUG > 3
    std::cout <<"dSdE ="<< dSdE;
#endif

    // Integrate the material stiffness matrix
    SymmTensor dEpsAi(nsd), dEpsBj(nsd);
    Matrix& Km = elMat.A[eKm-1];
    size_t a, b;
    unsigned short int i, j, k, l;
    for (b = 1; b <= Km.cols(); b++)
    {
      getDepsDu(dEpsBj,b);
      Matrix Ctmp(nsd,nsd);
      for (i = 1; i <= nsd; i++)
        for (j = 1; j <= nsd; j++)
          for (k = 1; k <= nsd; k++)
            for (l = 1; l <= nsd; l++)
              Ctmp(i,j) += dSdE(i,j,k,l)*dEpsBj(k,l);

      for (a = 1; a <= Km.rows(); a++)
      {
        getDepsDu(dEpsAi,a);
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


bool FractureElasticity::evalSol (Vector& s,
                                  const FiniteElement& fe, const Vec3& X,
                                  const std::vector<int>& MNPC) const
{
  // Extract element displacements
  Vectors eV(2);
  int ierr = 0;
  if (!primsol.empty() && !primsol.front().empty())
    ierr = utl::gather(MNPC,nsd,primsol.front(),eV.front());

  // Extract crack phase field vector for this element
  if (!myCVec.empty() && ierr == 0)
    ierr = utl::gather(MNPC,1,myCVec,eV.back());

  if (ierr > 0)
  {
    std::cerr <<" *** FractureElasticity::evalSol: Detected "<< ierr
              <<" node numbers out of range."<< std::endl;
    return false;
  }


  return this->evalSol2(s,eV,fe,X);
}


bool FractureElasticity::evalSol (Vector& s, const Vectors& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal, Vec3* pdir) const
{
  PROFILE3("FractureEl::evalSol");

  if (eV.size() < 2)
  {
    std::cerr <<" *** FractureElasticity::evalSol: Missing solution vector."
              << std::endl;
    return false;
  }
  else if (!eV.front().empty() && eV.front().size() != fe.dNdX.rows()*nsd)
  {
    std::cerr <<" *** FractureElasticity::evalSol: Invalid displacement vector."
              <<"\n     size(eV) = "<< eV.front().size() <<"   size(dNdX) = "
              << fe.dNdX.rows() <<","<< fe.dNdX.cols() << std::endl;
    return false;
  }
  else if (!eV[1].empty() && eV[1].size() != fe.N.size())
  {
    std::cerr <<" *** FractureElasticity::evalSol: Invalid phase field vector."
              <<"\n     size(eC) = "<< eV[1].size() <<"   size(N) = "
              << fe.N.size() << std::endl;
    return false;
  }

  // Evaluate the symmetric strain tensor, eps
  Matrix Bmat;
  SymmTensor eps(nsd);
  if (!this->kinematics(eV.front(),fe.N,fe.dNdX,0.0,Bmat,eps,eps))
    return false;

  // Evaluate the material parameters at this point
  double lambda, mu;
  if (!material->evaluate(lambda,mu,fe,X))
    return false;

  // Evaluate the stress state at this point
  SymmTensor sigma(nsd);
  double Phi, Gc = this->getStressDegradation(fe.N,eV[1]);
  if (!this->evalStress(lambda,mu,Gc,eps,Phi,sigma))
    return false;

  Vec3 p;
  bool havePval = false;
  if (toLocal && wantPrincipalStress)
  {
    // Calculate principal stresses and associated direction vectors
    havePval = pdir ? sigma.principal(p,pdir,2) : sigma.principal(p);
    // Congruence transformation to local coordinate system at current point
    if (locSys) sigma.transform(locSys->getTmat(X));
  }

  s = sigma;
  s.insert(s.begin()+2,0.2*(sigma(1,1)+sigma(2,2)));

  if (toLocal)
    s.push_back(sigma.vonMises());

  if (havePval)
  {
    s.push_back(p.x);
    s.push_back(p.y);
    if (sigma.dim() == 3)
      s.push_back(p.z);
  }

  if (toLocal)
  {
    s.push_back(Phi);
    s.push_back(Gc);
  }

  return true;
}


bool FractureElasticity::evalStress (double lambda, double mu, double Gc,
                                     const SymmTensor& epsilon,
                                     double& Phi, SymmTensor& sigma) const
{
  return this->evalStress(lambda,mu,Gc,epsilon,Phi,sigma,nullptr);
}


size_t FractureElasticity::getNoFields (int fld) const
{
  return this->Elasticity::getNoFields(fld) + (fld == 2 ? 2 : 0);
}


std::string FractureElasticity::getField2Name (size_t i, const char* pfx) const
{
  if (i < this->Elasticity::getNoFields(2))
    return this->Elasticity::getField2Name(i,pfx);

  std::string name;
  if (i == this->Elasticity::getNoFields(2))
    name = "Tensile energy, Phi";
  else
    name = "Stress degradation, g(c)";

  if (pfx) name = std::string(pfx) + " " + name;
  return name;
}
