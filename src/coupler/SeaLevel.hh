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

#ifndef SEALEVEL_H
#define SEALEVEL_H

#include "pism/util/Component.hh"

namespace pism {

class Geometry;

namespace ocean {
namespace sea_level {

class SeaLevel : public Component {
public:
  // "modifier" constructor
  SeaLevel(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> input);
  // "model" constructor
  SeaLevel(IceGrid::ConstPtr g);

  virtual ~SeaLevel();

  void init(const Geometry &geometry);

  void update(const Geometry &geometry, double t, double dt);

  const IceModelVec2S& elevation() const;

protected:
  virtual void init_impl(const Geometry &geometry);

  virtual void update_impl(const Geometry &geometry, double t, double dt);

  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

protected:
  std::shared_ptr<SeaLevel> m_input_model;
  IceModelVec2S m_sea_level;
};

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism

#endif /* SEALEVEL_H */
