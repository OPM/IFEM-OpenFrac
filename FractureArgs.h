// $Id$
//==============================================================================
//!
//! \file FractureArgs.h
//!
//! \date Jan 11 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Preparsing of input files for the FractureDynamics application.
//!
//==============================================================================

#ifndef _FRACTURE_ARGS_H
#define _FRACTURE_ARGS_H

#include "SIMargsBase.h"


/*!
  \brief Class holding the command-line argument values.
*/

class FractureArgs : public SIMargsBase
{
public:
  char* inpfile; //!< The input file to parse
  //! Time integrator to use (0=linear quasi-static, no phase-field coupling,
  //! 1=linear Newmark, 2=linear generalized alpha, 3=nonlinear quasi-static,
  //! 4=nonlinear Hilber-Hughes-Taylor)
  char integrator;
  char coupling; //!< Coupling flag (0: none, 1: staggered, 2: semi-implicit)
  bool poroEl;   //!< If \e true, use the poroelastic solver
  bool expPhase; //!< If \e true, use an explicit phase field

  //! \brief Default constructor.
  FractureArgs();

  //! \brief Parses a command-line argument.
  virtual bool parseArg(const char* argv);
  //! \brief Pre-parses an input file for command-line argument values.
  void parseFile(const char* argv, int& iarg);

protected:
  //! \brief Parse an element from the input file
  virtual bool parse(const TiXmlElement* elem);
};

#endif
