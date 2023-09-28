/* Copyright (C) 2023 PISM Authors
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

#include "pism/util/Component.hh"
#include <memory>
#include <vector>

namespace pism {

namespace array {
class Array3D;
class Array3DCollection;
class Scalar;
} // namespace array

namespace stressbalance {
class StressBalance;
}

/*!
 * The isochrone tracing scheme of [@ref Born2016] and [@ref Born2021].
 */
class Isochrones : public Component {
public:
  Isochrones(std::shared_ptr<const Grid> grid,
             std::shared_ptr<const stressbalance::StressBalance> stress_balance);
  virtual ~Isochrones() = default;

  void bootstrap(const array::Scalar &ice_thickness);

  void restart(const File &input_file, int record);

  void update(double t, double dt,
              const array::Array3D &u,
              const array::Array3D &v,
              const array::Scalar &ice_thickness,
              const array::Scalar &top_surface_mass_balance,
              const array::Scalar &bottom_surface_mass_balance);

  const array::Array3DCollection& layer_thicknesses() const;

private:
  MaxTimestep max_timestep_impl(double t) const;

  MaxTimestep max_timestep_cfl(double t) const;
  MaxTimestep max_timestep_deposition_times(double t) const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  double top_layer_deposition_time() const;

  DiagnosticList diagnostics_impl() const;

  void allocate(const std::vector<double> &levels);

  //! isochronal layer thicknesses
  std::shared_ptr<array::Array3DCollection> m_layer_thickness;

  //! temporary storage needed for time stepping
  std::shared_ptr<array::Array3DCollection> m_tmp;

  //! The index of the topmost isochronal layer.
  size_t m_top_layer_index;

  std::vector<double> m_deposition_times;

  std::shared_ptr<const stressbalance::StressBalance> m_stress_balance;
};

} // end of namespace pism
