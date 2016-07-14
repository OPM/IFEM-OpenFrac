// $Id$
//==============================================================================
//!
//! \file CubicMinimum.h
//!
//! \date Jul 13 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class finding the minimum of a cubic hermite interpolant.
//!
//==============================================================================

#include "CubicMinimum.h"
#include <GoTools/geometry/HermiteInterpolator.h>
#include <GoTools/geometry/SplineCurve.h>


bool CubicMinimum::Find (double& alpha,
                         const std::vector<double>& params,
                         const std::vector<double>& vals,
                         const std::vector<double>& tgts)
{
  if (vals.empty() || tgts.size() != vals.size())
    return false;

  // insert data in GoTools structure
  std::vector<Go::Point> samples;
  samples.reserve(2*vals.size());

  size_t i;
  Go::Point val(1);
  for (i = 0; i < vals.size(); i++) {
    val[0] = vals[i];
    samples.push_back(val);
    val[0] = tgts[i];
    samples.push_back(val);
  }

  // interpolate
  Go::HermiteInterpolator interp;
  std::vector<double> coefs;
  interp.interpolate(samples, params, coefs);

  // create curve
  Go::SplineCurve crv(interp.basis(), coefs.begin(), 1);

  // grab tangent curve
  Go::SplineCurve* dcrv = crv.derivCurve(1);

  // for each knot-span find point closest to zero and mark as an extremum,
  // if close enough
  val[0] = 0.0;
  std::vector<double> extrema;
  for (i = 0; i+1 < params.size(); i++) {
    Go::Point loc_pt;
    double loc_alpha, loc_dist;
    dcrv->closestPoint(val, params[i], params[i+1],
                       loc_alpha, loc_pt, loc_dist);
    if (loc_dist < 1.0e-5)
      extrema.push_back(loc_alpha);
  }
  delete dcrv;

  // no extrema found
  if (extrema.empty())
    return false;

  // choose acceptable solution among values
  dcrv = crv.derivCurve(2);
  alpha = 0.0;
  double minVal = 1.0e99;
  for (const double& it : extrema) {
    Go::Point pt_dd, pt;
    crv.point(pt, it);
    dcrv->point(pt_dd, it);
    if (pt_dd[0] >= 0.0 && pt[0] < minVal) {
      alpha = it;
      minVal = pt[0];
    }
  }

#ifdef SP_DEBUG
  std::cout <<"CubicMinimum::Find: alpha="<< alpha
            <<" minVal="<< minVal << std::endl;
#endif
  delete dcrv;
  return true;
}
