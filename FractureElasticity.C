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
#include "TimeDomain.h"
#include "MaterialBase.h"
#include "ElmMats.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Tensor4.h"
#include "Tensor.h"
#include "Profiler.h"
#include "IFEM.h"
#include "tinyxml.h"

#ifndef epsZ
//! \brief Zero tolerance for strains.
#define epsZ 1.0e-16
#endif


FractureElasticity::FractureElasticity (unsigned short int n)
  : Elasticity(n), mySol(primsol)
{
  crackP = nullptr;
  sigmaC = alpha = alpha1 = alpha0 = tSplit = 0.0;
  crpCut = zeta = 1.0;
  this->registerVector("phasefield",&myCVec);
  eC = 1; // Assuming second vector is phase field
}


FractureElasticity::FractureElasticity (IntegrandBase* parent,
                                        unsigned short int n)
  : Elasticity(n), mySol(parent->getSolutions())
{
  crackP = nullptr;
  sigmaC = alpha = alpha1 = alpha0 = tSplit = 0.0;
  crpCut = zeta = 1.0;
  parent->registerVector("phasefield",&myCVec);
  // Assuming second vector is pressure, third vector is pressure velocity
  eC = 3; // and fourth vector is the phase field
}


bool FractureElasticity::parse (const TiXmlElement* elem)
{
  const char* value = utl::getValue(elem,"stabilization");
  if (value)
  {
    alpha0 = alpha1 = atof(value);
    utl::getAttribute(elem,"initial",alpha0);
  }
  else if ((value = utl::getValue(elem,"crackpressure")))
  {
    std::string type;
    utl::getAttribute(elem,"type",type);
    utl::getAttribute(elem,"cutoff",crpCut);
    IFEM::cout <<"\tCrack pressure";
    crackP = utl::parseRealFunc(value,type);
    IFEM::cout << std::endl;
  }
  else if ((value = utl::getValue(elem,"critical_stress")))
  {
    sigmaC = atof(value);
    utl::getAttribute(elem,"slope",zeta);
  }
  else if ((value = utl::getValue(elem,"noEnergySplit")))
    tSplit = atof(value);
  else if (!strcasecmp(elem->Value(),"noEnergySplit"))
    tSplit = -1.0;
  else
    return this->Elasticity::parse(elem);

  return true;
}


void FractureElasticity::printLog () const
{
  this->Elasticity::printLog();

  if (crackP && crpCut < 1.0)
    IFEM::cout <<"\tCrack pressure cut-off value: "<< crpCut << std::endl;

  if (intPrm[3] > 0.0)
    IFEM::cout <<"\tNo geometric stiffness."<< std::endl;

  if (alpha1 != 0.0)
  {
    IFEM::cout <<"\tStabilization parameter: "<< alpha1;
    if (alpha0 > alpha1) IFEM::cout <<" ("<< alpha0 <<")";
    IFEM::cout << std::endl;
  }

  if (sigmaC > 0.0)
    IFEM::cout <<"\tCritical stress: "<< sigmaC
               <<" slope parameter: "<< zeta << std::endl;

  if (tSplit < 0.0)
    IFEM::cout <<"\tIsotropic degrading of";
  else if (tSplit > 0.0)
    IFEM::cout <<"\tIsotropic degrading of strain energy density up to t="
               << tSplit <<",\n\tthereafter degrading of tensile";
  else
    IFEM::cout <<"\tDegrading of tensile";

  IFEM::cout <<" strain energy density."<< std::endl;
}


void FractureElasticity::setMode (SIM::SolutionMode mode)
{
  this->Elasticity::setMode(mode);
  if (eC <= 1) return; // no parent integrand

  eKg = 0; // No geometric stiffness (assuming linear behaviour)
  eM = eS = 0; // Inertia and external forces are calculated by parent integrand
  if (eKm) eKm = 2; // Index for stiffness matrix in parent integrand
  if (iS) iS = 2; // Index for internal force vector in parent integrand
  eC = mode == SIM::DYNAMIC ? 5 : 3; // include velocity & acceleration vectors
}


void FractureElasticity::initIntegration (size_t nGp, size_t)
{
  // Initialize internal tensile energy buffer
  myPhi.resize(nGp);
}


void FractureElasticity::initIntegration (const TimeDomain& prm,
                                          const Vector&, bool)
{
  alpha = prm.first && alpha0 > alpha1 ? alpha0 : alpha1;
}


bool FractureElasticity::initElement (const std::vector<int>& MNPC,
                                      LocalIntegral& elmInt)
{
  if (mySol.empty())
  {
    std::cerr <<" *** FractureElasticity::initElement:"
              <<" No primary solution vectors."<< std::endl;
    return false;
  }

  int ierr = 0;
  // Unless elmInt.vec is empty, assume it already contains the element-level
  // displacement, velocity and acceleration vectors at this point
  if (elmInt.vec.empty())
  {
    elmInt.vec.resize(mySol.size()+eC);

    // Extract displacement vector for this element
    if (!mySol.front().empty())
      ierr = utl::gather(MNPC,npv,mySol.front(),elmInt.vec.front());

    // Extract velocity and acceleration vectors for this element
    for (size_t i = 1; i < mySol.size(); i++)
      if (ierr == 0 && !mySol[i].empty())
        ierr = utl::gather(MNPC,npv,mySol[i],elmInt.vec[eC+i]);
  }

  // Extract crack phase field vector for this element
  if (ierr == 0 && !myCVec.empty())
    ierr = utl::gather(MNPC,1,myCVec,elmInt.vec[eC]);

#if INT_DEBUG > 2
  for (size_t i = 0; i < elmInt.vec.size(); i++)
    std::cout <<"Element solution vector "<< i+1 << elmInt.vec[i];
#endif

  if (ierr == 0) return true;

  std::cerr <<" *** FractureElasticity::initElement: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


LocalIntegral* FractureElasticity::getLocalIntegral (size_t nen, size_t,
                                                     bool neumann) const
{
  LocalIntegral* li = this->Elasticity::getLocalIntegral(nen,0,neumann);
  if (m_mode >= SIM::RHS_ONLY && !neumann)
    static_cast<ElmMats*>(li)->c.resize(1); // Total strain energy
  return li;
}


double FractureElasticity::MieheCrit56 (const Vec3& eps,
                                        double lambda, double mu) const
{
  double D = 0.0;
  for (unsigned short int a = 0; a < nsd; a++)
    if (eps[a] > 0.0)
      D += eps[a]*eps[a];

  if (D == 0.0)
    return D;

  double E = mu+mu + lambda*mu/(lambda+mu); // Youngs modulus

  D *= E*E/(sigmaC*sigmaC);
  return D > 1.0 ? zeta*(D-1.0) : 0.0;
}


bool FractureElasticity::evalStress (double lambda, double mu, double Gc,
                                     const SymmTensor& epsilon,
                                     double* Phi, SymmTensor& sigma) const
{
  return this->evalStress(lambda,mu,Gc,epsilon,Phi,sigma,nullptr);
}


bool FractureElasticity::evalStress (double lambda, double mu, double Gc,
                                     const SymmTensor& epsilon, double* Phi,
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
    sigma = Phi[0] = Phi[1] = Phi[2] = 0.0;
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

  // Evaluate the stress tensor
  sigma = C0*trEps;
  sigma += 2.0*mu*(Gc*ePos + eNeg);

  // Evaluate the tensile energy
  Phi[0] = mu*(ePos*ePos).trace();
  if (trEps > 0.0) Phi[0] += 0.5*lambda*trEps*trEps;
  // Evaluate the compressive energy
  Phi[1] = mu*(eNeg*eNeg).trace();
  if (trEps < 0.0) Phi[1] += 0.5*lambda*trEps*trEps;
  // Evaluate the total strain energy
  Phi[2] = Gc*Phi[0] + Phi[1];

  if (sigmaC > 0.0) // Evaluate the Miehe crack driving function
    Phi[0] = this->MieheCrit56(eps,lambda,mu);

#if INT_DEBUG > 4
  std::cout <<"eps_p = "<< eps <<"\n";
  for (a = 0; a < nsd; a++)
    std::cout <<"M("<< 1+a <<") =\n"<< M[a];
  std::cout <<"ePos =\n"<< ePos <<"eNeg =\n"<< eNeg <<"sigma =\n"<< sigma
            <<"Phi = "<< Phi[0] <<" "<< Phi[1] <<" "<< Phi[2] << std::endl;
#endif

  if (!dSdE)
    return true;
  else if (eps[0] == eps[nsd-1])
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
                                                 const Vectors& eV) const
{
  // Evaluate the crack phase field function, c(X)
  double c = eV[eC].empty() ? 1.0 : N.dot(eV[eC]);
  if (c > 1.0) c = 1.0; // Ignore values larger than 1.0
  // Evaluate the stress degradation function, g(c), ignoring negative values
  return c > 0.0 ? (1.0-alpha)*c*c + alpha : alpha;
}


bool FractureElasticity::formCrackForce (Vector& ES, const Vectors& eV,
                                         const FiniteElement& fe,
                                         const Vec3& X) const
{
  if (!crackP)
    return true; // No crack pressure function, silently ignore

  if (crpCut < 1.0 && eV[eC].dot(fe.N) > crpCut)
    return true; // No crack pressure in the undamaged material

  Vector gradC; // Evaluate the phase field gradient, gradC = dNdX^t*eC
  if (!fe.dNdX.multiply(eV[eC],gradC,true))
    return false;

  // Integrate the load vector due to internal crack pressure
  double cpJW = (*crackP)(X) * fe.detJxW;
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += fe.N(a)*gradC(i)*cpJW;

  return true;
}


bool FractureElasticity::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe, const Vec3& X) const
{
  PROFILE3("FractureEl::evalInt");

  if (tSplit != 0.0)
  {
    std::cerr <<" *** FractureElasticity::evalInt: Isotropic degrading is not"
              <<" available with the tensorial formulation."<< std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  Matrix Bmat;
  Tensor4 dSdE(nsd);
  SymmTensor eps(nsd), sigma(nsd);
  bool lHaveStrains = false;
  double U = 0.0;

  if (eKm || eKg || iS || m_mode == SIM::RECOVERY)
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
    double Gc = this->getStressDegradation(fe.N,elmInt.vec);
#if INT_DEBUG > 3
    std::cout <<"lambda = "<< lambda <<" mu = "<< mu <<" G(c) = "<< Gc <<"\n";
    if (lHaveStrains) std::cout <<"eps =\n"<< eps;
#endif

    // Evaluate the stress state at this point
    double Phi[3];
    if (!this->evalStress(lambda,mu,Gc,eps,Phi,sigma, eKm ? &dSdE : nullptr))
      return false;

    if (m_mode != SIM::RHS_ONLY)
      myPhi[fe.iGP] = Phi[0];
    U = Phi[2];
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

  if (eS)
  {
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,fe.detJxW);
    // Integrate the load vector due to internal crack pressure
    if (!this->formCrackForce(elMat.b[eS-1],elMat.vec,fe,X))
      return false;
  }

  if (lHaveStrains && !elMat.c.empty())
    // Integrate the total strain energy
    elMat.c.front() += U*fe.detJxW;

  return true;
}


bool FractureElasticity::getElementSolution (Vectors& eV,
                                             const std::vector<int>& MNPC) const
{
  // Extract element displacements
  eV.resize(1+eC);
  int ierr = 0;
  if (!mySol.empty() && !mySol.front().empty())
    ierr = utl::gather(MNPC,nsd,mySol.front(),eV.front());

  // Extract crack phase field vector for this element
  if (!myCVec.empty() && ierr == 0)
    ierr = utl::gather(MNPC,1,myCVec,eV[eC]);

  if (ierr == 0)
    return true;

  std::cerr <<" *** FractureElasticity::getElementSolution: Detected "<< ierr
            <<" node numbers out of range."<< std::endl;
  return false;
}


bool FractureElasticity::evalSol (Vector& s,
                                  const FiniteElement& fe, const Vec3& X,
                                  const std::vector<int>& MNPC) const
{
  Vectors eV(1+eC);
  return this->getElementSolution(eV,MNPC) && this->evalSol2(s,eV,fe,X);
}


bool FractureElasticity::evalSol (Vector& s, const Vectors& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal, Vec3* pdir) const
{
  PROFILE3("FractureEl::evalSol");

  if (eV.size() <= eC)
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
  else if (!eV[eC].empty() && eV[eC].size() != fe.N.size())
  {
    std::cerr <<" *** FractureElasticity::evalSol: Invalid phase field vector."
              <<"\n     size(eC) = "<< eV[eC].size() <<"   size(N) = "
              << fe.N.size() << std::endl;
    return false;
  }

  // Evaluate the symmetric strain tensor, eps
  Matrix Bmat;
  SymmTensor eps(nsd);
  if (!this->kinematics(eV.front(),fe.N,fe.dNdX,0.0,Bmat,eps,eps))
    return false;
  else if (!eps.isZero(1.0e-16))
    for (unsigned short int i = 1; i <= nsd; i++)
      for (unsigned short int j = i+1; j <= nsd; j++)
        eps(i,j) *= 0.5; // Using tensor formulation

  // Evaluate the material parameters at this point
  double lambda, mu;
  if (!material->evaluate(lambda,mu,fe,X))
    return false;

  // Evaluate the stress state at this point
  SymmTensor sigma(nsd);
  double Phi[3];
  double Gc = this->getStressDegradation(fe.N,eV);
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
  if (nsd == 2)
  {
    // Insert the sigma_zz component for 2D plane strain
    double nu = 0.5*lambda/(lambda+mu);
    s.insert(s.begin()+2,nu*(sigma(1,1)+sigma(2,2)));
  }

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
    s.insert(s.end(),Phi,Phi+3);
    s.push_back(Gc);
  }

  return true;
}


double FractureElasticity::evalPhaseField (Vec3& gradD,
                                           const Vectors& eV,
                                           const Vector& N,
                                           const Matrix& dNdX) const
{
  if (eV.size() <= eC)
  {
    std::cerr <<" *** FractureElasticity::evalPhaseField:"
              <<" Missing phase field solution vector."<< std::endl;
    return -1.1;

  }
  else if (eV[eC].empty()) // No phase field ==> no crack yet
  {
    gradD = 0.0;
    return 0.0;
  }
  else if (eV[eC].size() != N.size())
  {
    std::cerr <<" *** FractureElasticity::evalPhaseField:"
              <<" Invalid phase field vector.\n     size(eC) = "
              << eV[eC].size() <<"   size(N) = "<< N.size() << std::endl;
    return -1.2;
  }

  // Invert the nodal phase field values, D = 1 - C,
  // since that is the convention used in the Miehe papers
  Vector eD(eV[eC].size()), tmp(nsd);
  for (size_t i = 0; i < eD.size(); i++)
    eD[i] = 1.0 - eV[eC][i];

  // Compute the phase field gradient, dD/dX
  if (dNdX.multiply(eD,tmp,true))
    gradD = tmp;
  else
    return -2.0;

  // Compute the phase field value, filtering out values outside [0.0,1.0]
  double d = eD.dot(N);
  return d > 1.0 ? 1.0 : (d < 0.0 ? 0.0 : d);
}


size_t FractureElasticity::getNoFields (int fld) const
{
  return this->Elasticity::getNoFields(fld) + (fld == 2 ? 4 : 0);
}


std::string FractureElasticity::getField2Name (size_t i, const char* pfx) const
{
  size_t nf = this->Elasticity::getNoFields(2);
  if (i < nf)
    return this->Elasticity::getField2Name(i,pfx);

  std::string name;
  if (i == nf)
    name = "Tensile energy density";
  else if (i == nf+1)
    name = "Compressive energy density";
  else if (i == nf+2)
    name = "Strain energy density";
  else
    name = "Stress degradation, g(c)";

  if (pfx) name = std::string(pfx) + " " + name;
  return name;
}
