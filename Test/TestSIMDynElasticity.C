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
#include "SIMElasticityWrap.h"

#include "NewmarkSIM.h"
#include "SIM2D.h"

#include <catch2/catch_test_macros.hpp>


TEST_CASE("TestSIMDynElasticity.Parse")
{
  SIMDynElasticity<SIM2D,NewmarkSIM,SIMElasticityWrap<SIM2D>> sim;
  REQUIRE(sim.read("Rectangle-p2.xinp"));
  REQUIRE(sim.preprocess());
}
