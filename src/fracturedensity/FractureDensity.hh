/* Copyright (C) 2019, 2021 PISM Authors
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
#ifndef FRACTUREDENSITY_H
#define FRACTUREDENSITY_H

#include "pism/util/IceGrid.hh"
#include "pism/util/Component.hh"
#include "pism/rheology/FlowLaw.hh"

namespace pism {

class IceModelVec2S;
class Geometry;

class FractureDensity : public Component {
public:
  FractureDensity(IceGrid::ConstPtr grid, std::shared_ptr<const rheology::FlowLaw> flow_law);
  virtual ~FractureDensity() = default;

  void restart(const File &input_file, int record);
  void bootstrap(const File &input_file);
  void initialize(const IceModelVec2S &density, const IceModelVec2S &age);
  void initialize();

  void update(double dt,
              const Geometry &geometry,
              const IceModelVec2V &velocity,
              const IceModelVec2S &hardness,
              const IceModelVec2S &inflow_boundary_mask);

  const IceModelVec2S& density() const;
  const IceModelVec2S& growth_rate() const;
  const IceModelVec2S& healing_rate() const;
  const IceModelVec2S& flow_enhancement() const;
  const IceModelVec2S& age() const;
  const IceModelVec2S& toughness() const;

private:

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  DiagnosticList diagnostics_impl() const;

  IceModelVec2S m_density;
  IceModelVec2S m_density_new;
  IceModelVec2S m_growth_rate;
  IceModelVec2S m_healing_rate;
  IceModelVec2S m_flow_enhancement;
  IceModelVec2S m_age;
  IceModelVec2S m_age_new;
  IceModelVec2S m_toughness;

  //! major and minor principal components of horizontal strain-rate tensor (temporary storage)
  IceModelVec3 m_strain_rates;

  //! components of horizontal stress tensor along axes and shear stress (temporary storage)
  IceModelVec3 m_deviatoric_stresses;

  //! Ghosted copy of the ice velocity
  IceModelVec2V m_velocity;

  std::shared_ptr<const rheology::FlowLaw> m_flow_law;
};

} // end of namespace pism

#endif /* FRACTUREDENSITY_H */
