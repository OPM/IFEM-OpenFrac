//==============================================================================
//!
//! \file TestSIMDynElasticity.C
//!
//! \date Nov 10 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for elasticity wrapper with crack phase field coupling.
//!
//==============================================================================

#include "SIMDynElasticity.h"
#include "SIM2D.h"
#include "SystemMatrix.h"
#include "QuasiStaticSIM.h"

#include "gtest/gtest.h"


/*!
  \brief Test wrapper class for QuasiStaticSIM.
*/

class TestQSSim : public QuasiStaticSIM
{
public:
  TestQSSim(SIMbase& sim) : QuasiStaticSIM(sim) {}
  virtual ~TestQSSim() {}

  int solveFirstIteration(TimeStep& tp)
  {
    std::cout <<"\n  step="<< tp.step <<"  time="<< tp.time.t << std::endl;

    tp.iter = 0;
    alpha = 1.0;
    model.setQuadratureRule(opt.nGauss[0],true);

    if (!model.updateDirichlet(tp.time.t,&solution.front()))
      return 1;
    if (!model.setMode(SIM::STATIC))
      return 2;
    if (!model.assembleSystem(tp.time,solution))
      return 3;
    if (!model.extractLoadVec(residual))
      return 4;
    if (!model.solveSystem(linsol,msgLevel-1))
      return 5;
    if (this->checkConvergence(tp) != SIM::OK)
      return 6;

    tp.iter++;
    if (!this->updateConfiguration(tp))
      return 7;
    if (!model.updateDirichlet())
      return 8;
    if (!model.assembleSystem(tp.time,solution))
      return 9;
    if (!model.solveSystem(linsol,msgLevel-1))
      return 10;

    return 0;
  }

  bool evalFDeriv(int ver, const TimeDomain& time, double& fDer)
  {
    version = ver;
    Vectors tmpSol;
    tmpSol.push_back(solution.front());
    if (version == 2) tmpSol.push_back(linsol);
    return this->evalEnergyFunctional(time,tmpSol,nullptr,&fDer);
  }
};


/*!
  \brief Test wrapper class for SIMDynElasticity.
*/

class TestSimDynEl : public SIMDynElasticity<SIM2D,TestQSSim>
{
public:
  TestSimDynEl() : SIMDynElasticity<SIM2D,TestQSSim>(true) {}
  virtual ~TestSimDynEl() {}

  int solveFirstIteration(TimeStep& tp)
  {
    return this->dSim.solveFirstIteration(tp);
  }

  bool evalFDeriv(int ver, const TimeStep& tp, double& fDer)
  {
    return this->dSim.evalFDeriv(ver,tp.time,fDer);
  }
};


TEST(TestSIMDynElasticity, Monolithic)
{
  IFEM::getOptions().discretization = ASM::LRSpline;
  Elastic::planeStrain = true;

  TestSimDynEl sim;
  ASSERT_TRUE(sim.read("Square-slit-p2.xinp"));
  ASSERT_TRUE(sim.preprocess());
  ASSERT_TRUE(sim.initSystem(SystemMatrix::SPARSE));

  TimeStep tp;
  tp.stopTime = 1.1;
  double fDer1, fDer2;
  ASSERT_TRUE(sim.init(tp));
  ASSERT_TRUE(sim.advanceStep(tp));
  ASSERT_TRUE(tp.increment());
  ASSERT_TRUE(sim.solveStep(tp));
  ASSERT_TRUE(sim.advanceStep(tp));
  ASSERT_TRUE(tp.increment());
  ASSERT_EQ(sim.solveFirstIteration(tp),0);

  FractureElasticNorm::extEnr = false;
  std::cout <<"\n  ** Evaluating f'(alpha) **"<< std::endl;
  ASSERT_TRUE(sim.evalFDeriv(1,tp,fDer1));
  std::cout <<"\tfDer1 = "<< fDer1 << std::endl;
  ASSERT_TRUE(sim.evalFDeriv(2,tp,fDer2));
  std::cout <<"\tfDer2 = "<< fDer2 << std::endl;
  ASSERT_NEAR(fDer1,fDer2,1.0e-8);
}
