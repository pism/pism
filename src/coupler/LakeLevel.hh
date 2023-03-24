/* Copyright (C) 2018, 2023 PISM Authors
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

#ifndef PISM_LAKE_LEVEL_H
#define PISM_LAKE_LEVEL_H

#include "pism/util/Component.hh"

namespace pism {

class Geometry;

namespace ocean {
namespace lake_level {

class LakeLevel : public Component {
public:
  // "modifier" constructor
  LakeLevel(IceGrid::ConstPtr g, std::shared_ptr<LakeLevel> input);
  // "model" constructor
  LakeLevel(IceGrid::ConstPtr g);

  virtual ~LakeLevel();

  void init(const Geometry &geometry);

  void update(const Geometry &geometry, double t, double dt);

  const IceModelVec2S& elevation() const;

  bool expandMargins() const;

protected:
  virtual void init_impl(const Geometry &geometry);

  virtual void update_impl(const Geometry &geometry, double t, double dt);

  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual bool expandMargins_impl() const;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

protected:
  std::shared_ptr<LakeLevel> m_input_model;
  IceModelVec2S m_lake_level;
  double m_fill_value;
};

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

#endif /* PISM_LAKE_LEVEL_H */
