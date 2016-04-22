// $Id$
//==============================================================================
//!
//! \file FractureArgs.C
//!
//! \date Jan 11 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Preparsing of input files for the FractureDynamics application.
//!
//==============================================================================

#include "FractureArgs.h"
#include "Utilities.h"
#include "tinyxml.h"


FractureArgs::FractureArgs () : SIMargsBase("fracturedynamics")
{
  inpfile = nullptr;
  integrator = coupling = 1;
  poroEl = false;
}


bool FractureArgs::parseArg (const char* argv)
{
  if (!strcmp(argv,"-nocrack"))
    coupling = 0;
  else if (!strcmp(argv,"-semiimplicit"))
    coupling = 2;
  else if (!strcmp(argv,"-lstatic"))
    integrator = 0;
  else if (!strcmp(argv,"-GA"))
    integrator = 2;
  else if (!strcmp(argv,"-qstatic"))
    integrator = 3;
  else if (!strcmp(argv,"-Miehe"))
    integrator = coupling = 3;
  else if (!strcmp(argv,"-HHT"))
    integrator = 4;
  else if (!strcmp(argv,"-oldHHT"))
    integrator = 5;
  else if (!strncmp(argv,"-poro",5))
    poroEl = true;
  else if (!strncmp(argv,"-explcr",7))
    expPhase = true;
  else if (!strncmp(argv,"-noadap",7))
    adap = false;
  else
    return this->SIMargsBase::parseArg(argv);

  return true;
}


void FractureArgs::parseFile (const char* argv, int& iarg)
{
  inpfile = const_cast<char*>(argv);
  if (strcasestr(inpfile,".xinp"))
    if (this->readXML(inpfile,false))
      iarg = 0; // start over and let command-line options override input file
}


bool FractureArgs::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"fracturedynamics")) {
    utl::getAttribute(elem,"timeintegrator",integrator);
    utl::getAttribute(elem,"coupling",coupling);
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"semiimplicit"))
        coupling = 2;
      else if (!strcasecmp(child->Value(),"generalizedalpha"))
        integrator = 2;
      else if (!strcasecmp(child->Value(),"quasistatic"))
        integrator = 3;
      else if (!strcasecmp(child->Value(),"miehe"))
        integrator = coupling = 3;
      else if (!strcasecmp(child->Value(),"hht"))
        integrator = 4;
      else if (!strcasecmp(child->Value(),"hilberhughestaylor"))
        integrator = 4;
      else if (!strcasecmp(child->Value(),"poroelastic"))
        poroEl = true;
  }

  return this->SIMargsBase::parse(elem);
}
