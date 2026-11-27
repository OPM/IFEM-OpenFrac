// $Id$
//==============================================================================
//!
//! \file SIMFractureDynamics.C
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for fracture-dynamic problems.
//!
//==============================================================================

#include "SIMFractureDynamics.h"

#include "Profiler.h"
#include "SIMadmin.h"


bool readModel (SIMadmin& model, const std::string& infile)
{
  PROFILE("Model input");
  return model.read(infile.c_str());
}
