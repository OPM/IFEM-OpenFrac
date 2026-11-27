// $Id$
//==============================================================================
//!
//! \file FractureLaunchImpl.h
//!
//! \date Nov 27 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Definitions of the launch templates declared in FractureLaunch.h.
//! \details Included only by the translation units that instantiate them,
//! FractureLaunch\<Integrator\>\<Dim\>D.C, so that the main program sees the
//! declarations only.
//!
//==============================================================================

#ifndef _FRACTURE_LAUNCH_IMPL_H
#define _FRACTURE_LAUNCH_IMPL_H

#include "FractureLaunch.h"


template<class Dim, class Integrator>
int runIntegrator (const FractureArgs& args)
{
  if constexpr (std::is_same_v<Integrator,NonLinSIM>)
    switch (args.coupling) {
    case 2:
      return runSimulator1<Dim,NonLinSIM,QstaticFracture>(args);
    case 3:
      return runSimulator1<Dim,NonLinSIM,MieheFracture>(args);
    default:
      return runSimulator1<Dim,NonLinSIM>(args);
    }
  else
    return runSimulator1<Dim,Integrator>(args);
}

#endif
