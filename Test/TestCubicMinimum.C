//==============================================================================
//!
//! \file TestCubicMinimum.C
//!
//! \date Jul 13 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for class finding the minimum of a cubic hermite interpolant.
//!
//==============================================================================

#include "CubicMinimum.h"

#include "gtest/gtest.h"


class CubicFunction {
public:
  static double value(double x) { return (x*x - x)*x; }
  static double tangent(double x) { return (3.0*x - 2.0)*x; }
};


class LinearFunction {
public:
  static double value(double x) { return x; }
  static double tangent(double) { return 1.0; }
};


class QuadraticFunction {
public:
  static double value(double x) { return (x+0.25)*(x+0.75); }
  static double tangent(double x) { return x+x+1.0; }
};


TEST(TestCubicMinimum, CubicFunction)
{
  double alpha;
  std::vector<double> params(10), vals(10), tgts(10);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = -1.0 + 2.0/9.0 * i;
    vals[i] = CubicFunction::value(params[i]);
    tgts[i] = CubicFunction::tangent(params[i]);
  }
  ASSERT_TRUE(CubicMinimum::Find(alpha, params, vals, tgts));
  ASSERT_FLOAT_EQ(alpha, 2.0/3.0);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = 1.0/9.0 * i;
    vals[i] = CubicFunction::value(params[i]);
    tgts[i] = CubicFunction::tangent(params[i]);
  }
  ASSERT_TRUE(CubicMinimum::Find(alpha, params, vals, tgts));
  ASSERT_FLOAT_EQ(alpha, 2.0/3.0);
}


TEST(TestCubicMinimum, QuadraticFunction)
{
  double alpha;
  std::vector<double> params(10), vals(10), tgts(10);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = -1.0 + 2.0/9.0 * i;
    vals[i] = QuadraticFunction::value(params[i]);
    tgts[i] = QuadraticFunction::tangent(params[i]);
  }
  ASSERT_TRUE(CubicMinimum::Find(alpha, params, vals, tgts));
  ASSERT_FLOAT_EQ(alpha, -0.5);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = 1.0/9.0 * i;
    vals[i] = QuadraticFunction::value(params[i]);
    tgts[i] = QuadraticFunction::tangent(params[i]);
  }
  ASSERT_FALSE(CubicMinimum::Find(alpha, params, vals, tgts));
}


TEST(TestCubicMinimum, LinearFunction)
{
  double alpha;
  std::vector<double> params(10), vals(10), tgts(10);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = -1.0 + 2.0/9.0 * i;
    vals[i] = LinearFunction::value(params[i]);
    tgts[i] = LinearFunction::tangent(params[i]);
  }
  ASSERT_FALSE(CubicMinimum::Find(alpha, params, vals, tgts));

  for (size_t i = 0; i < 10; ++i) {
    params[i] = 1.0/9.0 * i;
    vals[i] = LinearFunction::value(params[i]);
    tgts[i] = LinearFunction::tangent(params[i]);
  }
  ASSERT_FALSE(CubicMinimum::Find(alpha, params, vals, tgts));
}


TEST(TestCubicMinimum, DiscreteValues)
{
  double alpha;
  std::vector<double> params(10);
  for (size_t i = 0; i < 10; ++i)
    params[i] = -1.0 + 2.0/9.0 * i;

  std::vector<double> vals({ 0.668273, 0.668279, 0.668284, 0.66829, 0.668296,
			     0.668301, 0.668307, 0.668311, 0.668306, 0.6683 });

  std::vector<double> tgts({ 1.15625e-06, 1.00442e-06, 8.52527e-07, 7.00576e-07, 5.48567e-07,
			     3.96498e-07, 2.44365e-07, 8.6514e-08, -1.15779e-07, -3.17837e-07 });

  ASSERT_TRUE(CubicMinimum::Find(alpha, params, vals, tgts));
  ASSERT_NEAR(alpha, -1.0, 1.0e-8);
}
