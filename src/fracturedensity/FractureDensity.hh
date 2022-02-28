/* Copyright (C) 2019, 2021, 2022 PISM Authors
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
#include "pism/util/array/Scalar.hh"
#include "pism/util/IceModelVec2V.hh"
#include "pism/util/array/Array2D.hh"
#include "pism/stressbalance/StressBalance.hh"

namespace pism {

class Geometry;

class FractureDensity : public Component {
public:
  FractureDensity(IceGrid::ConstPtr grid, std::shared_ptr<const rheology::FlowLaw> flow_law);
  virtual ~FractureDensity() = default;

  void restart(const File &input_file, int record);
  void bootstrap(const File &input_file);
  void initialize(const array::Scalar &density, const array::Scalar &age);
  void initialize();

  void update(double dt,
              const Geometry &geometry,
              const IceModelVec2V &velocity,
              const array::Scalar &hardness,
              const array::Scalar &inflow_boundary_mask);

  const array::Scalar& density() const;
  const array::Scalar& growth_rate() const;
  const array::Scalar& healing_rate() const;
  const array::Scalar& flow_enhancement() const;
  const array::Scalar& age() const;
  const array::Scalar& toughness() const;

private:

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  DiagnosticList diagnostics_impl() const;

  array::Scalar1 m_density;
  array::Scalar m_density_new;
  array::Scalar m_growth_rate;
  array::Scalar m_healing_rate;
  array::Scalar m_flow_enhancement;
  array::Scalar1 m_age;
  array::Scalar m_age_new;
  array::Scalar m_toughness;

  //! major and minor principal components of horizontal strain-rate tensor (temporary storage)
  array::Array2D<stressbalance::PrincipalStrainRates> m_strain_rates;

  //! components of horizontal stress tensor along axes and shear stress (temporary storage)
  array::Array3D m_deviatoric_stresses;

  //! Ghosted copy of the ice velocity
  IceModelVec2V m_velocity;

  std::shared_ptr<const rheology::FlowLaw> m_flow_law;
};

} // end of namespace pism

#endif /* FRACTUREDENSITY_H */
