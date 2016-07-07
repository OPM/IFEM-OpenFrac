// $Id$
//==============================================================================
//!
//! \file SIMExplPhaseField.C
//!
//! \date May 27 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver representing an explicit phase-field.
//!
//==============================================================================

#include "SIMExplPhaseField.h"
#include "SIMoutput.h"
#include "TimeStep.h"
#include "DataExporter.h"
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


SIMExplPhaseField::SIMExplPhaseField (SIMoutput* gridOwner)
{
  myHeading = "Explicit phase field";
  myOwner   = gridOwner;
  phaseFunc = nullptr;
  myStep    = 0;
}


void SIMExplPhaseField::registerFields (DataExporter& exporter)
{
  exporter.registerField("c","phase field",DataExporter::SIM,
                         DataExporter::PRIMARY);
  exporter.setFieldValue("c",this,&phaseField);
}


bool SIMExplPhaseField::init (const TimeStep&)
{
  this->registerField("phasefield",phaseField);
  return true;
}


bool SIMExplPhaseField::parse (const TiXmlElement* elem)
{
  const char* value = utl::getValue(elem,"phasefield");
  if (!value)
    return this->SIMbase::parse(elem);

  std::string type;
  utl::getAttribute(elem,"type",type);
  IFEM::cout <<"\tPhase-field function";
  phaseFunc = utl::parseRealFunc(value,type);
  IFEM::cout << std::endl;

  return true;
}


bool SIMExplPhaseField::solveStep (TimeStep& tp, bool)
{
  if (!phaseFunc)
  {
    std::cerr <<" *** SIMExplPhaseField::solveStep: No phase field function."
              << std::endl;
    return false;
  }

  phaseField.resize(myOwner->getNoNodes(1));
  return myOwner->project(phaseField,phaseFunc,1,0,1,
                          SIMoptions::GLOBAL,tp.time.t);
}


bool SIMExplPhaseField::saveStep (const TimeStep& tp, int& nBlock)
{
  if (tp.step%opt.saveInc == 0 && opt.format >= 0)
    if (myOwner->writeGlvS1(phaseField,++myStep,nBlock,tp.time.t,
                            "phase",6,1,true) < 0) return false;

  return true;
}
