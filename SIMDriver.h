// $Id$
//==============================================================================
//!
//! \file SIMDriver.h
//!
//! \date Nov 15 2015
//!
//! \author Knut Morten Okstad
//!
//! \brief Simulation driver class template.
//!
//==============================================================================

#ifndef _SIM_DRIVER_H
#define _SIM_DRIVER_H


/*!
  \brief Simulation driver template class.

  \details Only the parse method is reimplemented here to handle that the
  time stepping parameters may be located within the specified context, and
  for obtaining the file name for global energy output from child simulator.
*/

template<class T, template<class S> class Solver>
class SIMDriver : public Solver<T>
{
public:
  //! \brief The constructor initializes the reference to the actual solver.
  SIMDriver(T& s, const char* c = nullptr) : Solver<T>(s), context(c) {}
  //! \brief Empty destructor.
  virtual ~SIMDriver() {}

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),context))
    {
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        this->SIMSolver<T>::parse(child);
    }
    else if (!strcasecmp(elem->Value(),"postprocessing"))
    {
      const TiXmlElement* child = elem->FirstChildElement("energyfile");
      if (child && child->FirstChild())
        this->S1.setEnergyFile(child->FirstChild()->Value());
    }

    return this->Solver<T>::parse(elem);
  }

private:
  const char* context; //!< XML-tag to search for time-stepping input within
};

#endif
