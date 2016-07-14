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
#include <cmath>


class CubicFunction {
public:
  static double value(double alpha)
  {
    return pow(alpha, 3.0) - pow(alpha, 2.0);
  }

  static double tangent(double alpha)
  {
    return 3.0*pow(alpha,2.0) - 2.0*alpha;
  }
};


class LinearFunction {
public:
  static double value(double alpha)
  {
    return alpha;
  }

  static double tangent(double alpha)
  {
    return 1.0;
  }
};


class QuadraticFunction {
public:
  static double value(double alpha)
  {
    return (alpha+0.25)*(alpha+0.75);
  }

  static double tangent(double alpha)
  {
    return 2*alpha + 1;
  }
};


TEST(TestCubicMinimum, CubicFunction)
{
  double alpha;
  std::vector<double> params(10);
  std::vector<double> vals(10);
  std::vector<double> tgts(10);
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
  std::vector<double> params(10);
  std::vector<double> vals(10);
  std::vector<double> tgts(10);
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
  std::vector<double> params(10);
  std::vector<double> vals(10);
  std::vector<double> tgts(10);
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
