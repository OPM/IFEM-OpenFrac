// $Id$
//==============================================================================
//!
//! \file SIMCH.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for Cahn-Hilliard problems.
//!
//==============================================================================

#ifndef _SIM_CH_H
#define _SIM_CH_H

#include "CahnHilliard.h"
#include "AnaSol.h"
#include "ASMstruct.h"
#include "Functions.h"
#include "InitialConditionHandler.h"
#include "Property.h"
#include "SIMoutput.h"
#include "SIMSolver.h"
#include "IFEM.h"
#include "TimeStep.h"
#include "Profiler.h"
#include "Utilities.h"
#include "DataExporter.h"
#include "tinyxml.h"

#include <memory>


/*!
  \brief Driver class for an Cahn-Hilliard simulator.
*/

template<class Dim, class Integrand=CahnHilliard> class SIMCH : public Dim
{
public:
  typedef bool SetupProps;

  //! \brief Default constructor.
  SIMCH() : Dim(1), ch(Dim::dimension), filter(false)
  {
    Dim::myProblem = &ch;
    this->myHeading = "Cahn-Hilliard solver";
  }

  //! \brief Destructor
  virtual ~SIMCH()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

  //! \brief Parses a data section from an input stream (deprecated).
  virtual bool parse(char*, std::istream&) { return false; }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"cahnhilliard"))
      return this->Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement()) {
      if (strcmp(child->Value(),"smearing") == 0) {
        const char* value = utl::getValue(child,"smearing");
        if (value)
          ch.setSmearFactor(atof(value));
      } else if (strcmp(child->Value(),"maxcrack") == 0) {
        const char* value = utl::getValue(child,"maxcrack");
        if (value)
          ch.setMaxCrack(atof(value));
      } else if (strcmp(child->Value(),"filter_values") == 0) {
        const static std::map<std::string,SIMoptions::ProjectionMethod> typemap =
                {{"global", SIMoptions::GLOBAL},
                 {"dgl2",   SIMoptions::DGL2},
                 {"cgl2",   SIMoptions::CGL2},
                 {"scr",    SIMoptions::SCR},
                 {"vdsa",   SIMoptions::VDSA},
                 {"quasi",  SIMoptions::QUASI},
                 {"lsq",    SIMoptions::LEASTSQ}};
        std::string type="global";
        utl::getAttribute(child,"type",type);
        auto it = typemap.find(type);
        if (it == typemap.end()) {
          it = typemap.begin();
          std::cerr << "Unknown projection method. Defaulting to global." << std::endl;
        }

        method = it->second;
        IFEM::cout << "\tFiltering phase field using " << it->first << " projection." << std::endl;
        filter = true;
      } else if (strcmp(child->Value(),"initial_crack") == 0) {
        std::string type;
        utl::getAttribute(child,"type",type);
        const char* value = utl::getValue(child,"initial_crack");
        if (value) {
          IFEM::cout << "\tInitial crack function";
          initial_crack.reset(utl::parseRealFunc(value,type));
        }
      } else if (!strcasecmp(child->Value(),"material")) {
        int code = this->parseMaterialSet(child,mVec.size());
        std::cout <<"\tMaterial code "<< code <<":" << std::endl;
        mVec.push_back(std::unique_ptr<CHMaterial>(new CHMaterial));
        mVec.back()->parse(child);
        mVec.back()->printLog();
        ch.setMaterial(mVec.back().get()); // in order to inject common material/for single patch models
      } else
        this->Dim::parse(child);
    }

    return true;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "CahnHilliard"; }

  //! \brief Initialize problem.
  bool init(const TimeStep& tp)
  {
    IFEM::cout << std::endl << "== Cahn-Hilliard problem setup ==" << std::endl;
    ch.printLog();
    if (this->hasIC("phasefield"))
      std::cout << "Initial phase field specified." << std::endl;
    else if (initial_crack) {
      IFEM::cout << "Initial crack specified as a function." << std::endl;
      ch.setInitialCrackFunction(initial_crack.get());
    }

    this->registerField("phasefield", phasefield);
    return true;
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  virtual bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Dummy function.
  virtual bool advanceStep(TimeStep& tp) { return true; }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMCH::solveStep");

    if (!this->assembleSystem())
      return false;

    if (!this->solveSystem(phasefield, Dim::msgLevel-1,"phasefield  "))
      return false;

    ch.setInitialCrackFunction(nullptr);

    return !filter || postSolve(tp,true);
  }

  //! \brief Make sure phase field value is between 0 and 1.
  bool postSolve(const TimeStep& tp,bool)
  {
    Matrix tmp;
    if (!this->project(tmp,phasefield,SIMoptions::DGL2))
      return false;

    phasefield = tmp.getRow(1);
    return true;
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMCH::saveStep");

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;

    // Write solution fields
    bool result = this->writeGlvS(phasefield, iDump, nBlock,
                                  tp.time.t, true, "phasefield", 89);

    return this->writeGlvStep(iDump, tp.time.t) && result;
  }

  //! \brief Register fields for data output.
  void registerFields(DataExporter& exporter, const std::string& prefix="")
  {
    exporter.registerField("c","phase field",DataExporter::SIM,
                           DataExporter::PRIMARY|DataExporter::RESTART,
                           prefix);
    exporter.setFieldValue("c", this, &phasefield);
  }

  //! \brief Sets initial conditions.
  void setInitialConditions() { SIM::setInitialConditions(*this); }

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd)
  {
    if (propInd >= mVec.size())
      propInd = mVec.size()-1;

    ch.setMaterial(mVec[propInd].get());
    return true;
  }

  //! \brief Set tensile energy vector from structure problem.
  void setTensileEnergy(const Vector* tens) { ch.setTensileEnergy(tens); }
private:
  Integrand ch; //!< Problem definition.
  std::unique_ptr<RealFunc> initial_crack; //!< Function describing an initial crack.
  std::vector<std::unique_ptr<CHMaterial>> mVec;  //!< Material data

  Vector phasefield; //!< Current phase field solution.
  bool filter; //!< Whether or not to filter results.
  SIMoptions::ProjectionMethod method; //!< Projection method used for filter.
};


//! \brief Partial specialization for configurator
template<class Dim, class Integrand>
struct SolverConfigurator< SIMCH<Dim,Integrand> > {
  int setup(SIMCH<Dim,Integrand>& ch,
            const typename SIMCH<Dim,Integrand>::SetupProps& props, char* infile)
  {
    utl::profiler->start("Model input");

    // Reset the global element and node numbers
    ASMstruct::resetNumbering();
    if (!ch.read(infile))
      return 2;

    utl::profiler->stop("Model input");

    // Preprocess the model and establish data structures for the algebraic system
    if (!ch.preprocess())
      return 3;

    // Initialize the linear solvers
    ch.setMode(SIM::DYNAMIC);
    ch.initSystem(ch.opt.solver,1,1,false);
    ch.setQuadratureRule(ch.opt.nGauss[0],true);

    // Time-step loop
    ch.init(TimeStep());
    ch.setInitialConditions();

    return 0;
  }
};

#endif
