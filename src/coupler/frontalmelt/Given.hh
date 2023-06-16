// Copyright (C) 2019 PISM Authors
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

#ifndef _PFM_GIVEN_H_
#define _PFM_GIVEN_H_

#include "pism/coupler/FrontalMelt.hh"

#include "pism/util/array/Forcing.hh"

namespace pism {
namespace frontalmelt {

class Given : public FrontalMelt
{
public:
  Given(std::shared_ptr<const Grid> g);
  virtual ~Given() = default;

private:
  void init_impl(const Geometry &geometry);

  void update_impl(const FrontalMeltInputs &inputs, double t, double dt);

  MaxTimestep max_timestep_impl(double t) const;

  const array::Scalar& frontal_melt_rate_impl() const;

  std::shared_ptr<array::Forcing> m_frontal_melt_rate;
};

} // end of namespace frontalmelt
} // end of namespace pism
#endif /* _PFM_GIVEN_H_ */
