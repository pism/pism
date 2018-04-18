// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _DELTA_SL_2D_
#define _DELTA_SL_2D_

#include "pism/coupler/SeaLevel.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {

class ScalarForcing;

namespace ocean {

namespace sea_level {

/*!
 * 2D sea level forcing
 */
class Delta_SL_2D : public SeaLevel {
public:
  Delta_SL_2D(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> in);
  virtual ~Delta_SL_2D();

private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  IceModelVec2T::Ptr m_forcing;
};

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism

#endif /* _DELTA_SL_2D_ */
