// $Id$
//==============================================================================
//!
//! \file FractureLaunchHHTSIM3D.C
//!
//! \date Nov 27 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Instantiates the launch of the HHTSIM time integrator for 3D
//! models.
//! \details One of the translation units the launch is split into, which
//! compile in parallel.
//!
//==============================================================================

#include "FractureLaunchImpl.h"

//! \cond
template int runIntegrator<SIM3D,HHTSIM>(const FractureArgs&);
//! \endcond
