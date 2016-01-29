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

#include "gtest/gtest.h"

TEST(TestSIMDynElasticity, Parse)
{
  SIMDynElasticity<SIM2D> sim;
  EXPECT_TRUE(sim.read("Rectangle-p2.xinp"));
  EXPECT_TRUE(sim.preprocess());
}
