//==============================================================================
//!
//! \file TestFractureElasticity.C
//!
//! \date Dec 15 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for the fracture elasticity integrands.
//!
//==============================================================================

#include "FractureElasticityVoigt.h"
#include "Tensor4.h"
#include "Tensor.h"
#include <iostream>

#include "gtest/gtest.h"


/*!
  \brief Wrapper class, needed because evalStress is (and should be) protected.
*/

class FracEl : public FractureElasticityVoigt
{
public:
  FracEl(unsigned short int n) : FractureElasticityVoigt(n) {}
  virtual ~FracEl() {}
  bool calcStress(double lambda, double mu, double Gc, const SymmTensor& eps,
                  double& Phi, SymmTensor& sigma, Matrix& dSdE) const
  { return evalStress(lambda,mu,Gc,eps,Phi,sigma,&dSdE); }
  bool calcStress(double lambda, double mu, double Gc, const SymmTensor& eps,
                  double& Phi, SymmTensor& sigma, Tensor4& dSdE) const
  { return FractureElasticity::evalStress(lambda,mu,Gc,eps,Phi,sigma,&dSdE); }
};


TEST(TestFractureElasticity, evalStress)
{
  size_t nsd = 2;
  size_t nstrc = (nsd+1)*nsd/2;
  FracEl frel(nsd);
  double lambda = 100.0, mu = 150.0, Gc = 1.0;
  SymmTensor eps(nsd), sigma(nsd);
  Matrix Cmat(nstrc,nstrc);
  Tensor4 dSdE(nsd);
  double Phi = 0.0;

  // Check that the Tensor- and matrix formulations give
  // equivalent constitutive matrix in the unloaded case
  EXPECT_TRUE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,Cmat));
  EXPECT_TRUE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,dSdE));
  std::cout <<"Cmat:"<< Cmat;
  std::cout <<"dSdE:"<< dSdE;
  EXPECT_FLOAT_EQ(Cmat(1,1),dSdE(1,1,1,1));
  EXPECT_FLOAT_EQ(Cmat(1,2),dSdE(1,1,2,2));
  EXPECT_FLOAT_EQ(Cmat(1,3),0.0);
  EXPECT_FLOAT_EQ(Cmat(2,2),dSdE(2,2,2,2));
  EXPECT_FLOAT_EQ(Cmat(2,3),0.0);
  EXPECT_FLOAT_EQ(Cmat(3,3),dSdE(1,2,1,2));

  // Check that the Tensor- and matrix formulations give
  // equivalent constitutive matrix in case of hydrostatic strain
  eps = 1.0;
  Cmat.fill(0.0);
  EXPECT_TRUE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,Cmat));
  EXPECT_TRUE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,dSdE));
  std::cout <<"Cmat:"<< Cmat;
  std::cout <<"dSdE:"<< dSdE;
  EXPECT_FLOAT_EQ(Cmat(1,1),dSdE(1,1,1,1));
  EXPECT_FLOAT_EQ(Cmat(1,2),dSdE(1,1,2,2));
  EXPECT_FLOAT_EQ(Cmat(1,3),0.0);
  EXPECT_FLOAT_EQ(Cmat(2,2),dSdE(2,2,2,2));
  EXPECT_FLOAT_EQ(Cmat(2,3),0.0);
  EXPECT_FLOAT_EQ(Cmat(3,3),dSdE(1,2,1,2));

  // Check that the Tensor- and matrix formulations give
  // equivalent constitutive matrix in case of uniaxial strain
  eps(2,2) = 0.0;
  Cmat.fill(0.0);
  EXPECT_TRUE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,Cmat));
  EXPECT_TRUE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,dSdE));
  std::cout <<"Cmat:"<< Cmat;
  std::cout <<"dSdE:"<< dSdE;
  EXPECT_NEAR(Cmat(1,1),dSdE(1,1,1,1),1.0e-8);
  EXPECT_NEAR(Cmat(1,2),dSdE(1,1,2,2),1.0e-8);
  EXPECT_NEAR(Cmat(1,3),0.0,1.0e-8);
  EXPECT_NEAR(Cmat(2,2),dSdE(2,2,2,2),1.0e-8);
  EXPECT_NEAR(Cmat(2,3),0.0,1.0e-8);
  EXPECT_NEAR(Cmat(3,3),dSdE(1,2,1,2),1.0e-8);

  // Check that the Tensor- and matrix formulations give
  // equivalent constitutive matrix in case of arbitrary strain
  eps(2,2) = 2.0;
  eps(1,2) = 0.5;
  Cmat.fill(0.0);
  EXPECT_TRUE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,Cmat));
  EXPECT_TRUE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,dSdE));
  std::cout <<"Cmat:"<< Cmat;
  std::cout <<"dSdE:"<< dSdE;
  EXPECT_NEAR(Cmat(1,1),dSdE(1,1,1,1),1.0e-8);
  EXPECT_NEAR(Cmat(1,2),dSdE(1,1,2,2),1.0e-8);
  EXPECT_NEAR(Cmat(1,3),0.0,1.0e-8);
  EXPECT_NEAR(Cmat(2,2),dSdE(2,2,2,2),1.0e-8);
  EXPECT_NEAR(Cmat(2,3),0.0,1.0e-8);
  EXPECT_NEAR(Cmat(3,3),dSdE(1,2,1,2),1.0e-8);
}
