// $Id$
//==============================================================================
//!
//! \file FractureLaunchNonLinSIM2D.C
//!
//! \date Nov 27 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Instantiates the launch of the NonLinSIM time integrator for 2D
//! models.
//! \details One of the translation units the launch is split into, which
//! compile in parallel.
//!
//==============================================================================

#include "FractureLaunchImpl.h"

//! \cond
template int runIntegrator<SIM2D,NonLinSIM>(const FractureArgs&);
//! \endcond
