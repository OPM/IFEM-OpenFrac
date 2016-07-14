// $Id$
//==============================================================================
//!
//! \file CubicMinimum.h
//!
//! \date Jul 13 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Functions for minimization of cubic interpolants.
//!
//==============================================================================

#include <vector>


namespace CubicMinimum //! Line-search utilities
{
  //! \brief Finds the minimum of a cubic hermite interpolant.
  bool Find(double& alpha,
            const std::vector<double>& params,
            const std::vector<double>& vals,
            const std::vector<double>& tgts);
}
