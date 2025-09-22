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

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinRel;


/*!
  \brief Wrapper class, needed because evalStress is (and should be) protected.
*/

class FracEl : public FractureElasticityVoigt
{
public:
  FracEl(unsigned short int n) : FractureElasticityVoigt(n) {}
  virtual ~FracEl() {}
  bool calcStress(double lambda, double mu, double Gc, const SymmTensor& eps,
                  double* Phi, SymmTensor& sigma, Matrix& dSdE) const
  { return evalStress(lambda,mu,Gc,eps,Phi,&sigma,&dSdE); }
  bool calcStress(double lambda, double mu, double Gc, const SymmTensor& eps,
                  double* Phi, SymmTensor& sigma, Tensor4& dSdE) const
  { return FractureElasticity::evalStress(lambda,mu,Gc,eps,Phi,sigma,&dSdE); }
};


static void compare (const Matrix& Cmat, const Tensor4& dSdE, size_t nsd)
{
  std::cout <<"Cmat:"<< Cmat;
  std::cout <<"dSdE:"<< dSdE;

  size_t i, j;
  for (i = 1; i <= nsd; i++)
    for (j = 1; j <= nsd; j++)
      REQUIRE_THAT(Cmat(i,j), WithinRel(dSdE(i,i,j,j)));

  for (i = 1; i <= nsd; i++)
    for (j = 1; nsd+j < Cmat.cols(); j++)
      REQUIRE_THAT(Cmat(i,nsd+j), WithinRel(dSdE(i,i,j,j%nsd+1)));

  for (i = 1; nsd+i <= Cmat.rows(); i++)
    for (j = 1; nsd+j <= Cmat.cols(); j++)
      REQUIRE_THAT(Cmat(nsd+i,nsd+j), WithinRel(dSdE(i,i%nsd+1,j,j%nsd+1)));
}


static void checkStress (size_t nsd)
{
  size_t nstrc = (nsd+1)*nsd/2;
  FracEl frel(nsd);
  double lambda = 100.0, mu = 150.0, Gc = 0.1;
  SymmTensor eps(nsd), sigma(nsd);
  Matrix Cmat(nstrc,nstrc);
  Tensor4 dSdE(nsd);
  double Phi[4];

  // Check that the Tensor- and matrix formulations give
  // equivalent constitutive matrix in the unloaded case
  REQUIRE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,Cmat));
  REQUIRE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,dSdE));
  compare(Cmat,dSdE,nsd);

  // Check that the Tensor- and matrix formulations give
  // equivalent constitutive matrix in case of hydrostatic strain
  eps = 1.0;
  Cmat.fill(0.0);
  REQUIRE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,Cmat));
  REQUIRE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,dSdE));
  compare(Cmat,dSdE,nsd);

  // Check that the Tensor- and matrix formulations give
  // equivalent constitutive matrix in case of uniaxial strain
  eps(2,2) = 0.0;
  Cmat.fill(0.0);
  REQUIRE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,Cmat));
  REQUIRE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,dSdE));
  compare(Cmat,dSdE,nsd);

  // Check that the Tensor- and matrix formulations give
  // equivalent constitutive matrix in case of arbitrary strain
  eps(2,2) = -0.5;
  eps(1,2) = 0.3;
  Cmat.fill(0.0);
  REQUIRE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,Cmat));
  REQUIRE(frel.calcStress(lambda,mu,Gc,eps,Phi,sigma,dSdE));
  compare(Cmat,dSdE,nsd);
}


TEST_CASE("TestFractureElasticity.EvalStress")
{
  const int param = GENERATE(2,3);
  SECTION("NSD = " + std::to_string(param)) {
    checkStress(param);
  }
}
