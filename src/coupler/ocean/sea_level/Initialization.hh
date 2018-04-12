/* Copyright (C) 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _SEA_LEVEL_INITIALIZATION_H_
#define _SEA_LEVEL_INITIALIZATION_H_

#include "pism/coupler/SeaLevel.hh"

namespace pism {
namespace ocean {
namespace sea_level {

/*! Sea level forcing "modifier" that helps with initialization.
 *
 * This modifier saves the sea level as a part of the model state and re-loads it during
 * initialization so that it is available *before* the first time step in a re-started
 * run.
 */
class InitializationHelper : public SeaLevel {
public:
  InitializationHelper(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> in);

private:
  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  void update_impl(const Geometry &geometry, double t, double dt);
  void init_impl(const Geometry &geometry);

  const IceModelVec2S& sea_level_elevation_impl() const;
};

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism

#endif /* _SEA_LEVEL_INITIALIZATION_H_ */
