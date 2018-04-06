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

#ifndef _PS_DELTA_T_H_
#define _PS_DELTA_T_H_

#include <memory>               // std::unique_ptr

#include "pism/coupler/SurfaceModel.hh"

namespace pism {

class ScalarForcing;

namespace surface {

/** @brief Implements the scalar temperature offsets for the ice
 * surface temperature.
 *
 * Other fields are passed through without change.
 */
class Delta_T : public SurfaceModel
{
public:
  Delta_T(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> in);
  virtual ~Delta_T();
protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  virtual const IceModelVec2S& temperature_impl() const;

  std::unique_ptr<ScalarForcing> m_forcing;

  IceModelVec2S::Ptr m_temperature;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PS_DELTA_T_H_ */
