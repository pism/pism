// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _PFMCONSTANT_H_
#define _PFMCONSTANT_H_

#include "pism/coupler/FrontalMelt.hh"

namespace pism {
namespace frontalmelt {

//! @brief A class implementing a constant (in terms of the ocean inputs) frontal melt
//! model. Uses a configuration parameter for the frontal melt rate.
class Constant : public FrontalMelt {
public:
  Constant(IceGrid::ConstPtr g);
  virtual ~Constant() = default;

private:
  void init_impl(const Geometry &geometry);

  void update_impl(const FrontalMeltInputs &inputs, double t, double dt);

  const array::Scalar& frontal_melt_rate_impl() const;

  MaxTimestep max_timestep_impl(double t) const;

  array::Scalar::Ptr m_frontal_melt_rate;
};

} // end of namespace frontalmelt
} // end of namespace pism
#endif /* _PFMCONSTANT_H_ */
