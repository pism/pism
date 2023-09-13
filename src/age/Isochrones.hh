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
#include "pism/util/array/Array3D.hh"
#include <memory>

namespace pism {

class Isochrones : public Component {
public:
  Isochrones(std::shared_ptr<const Grid> grid);
  virtual ~Isochrones() = default;

  void init(const Geometry &geometry);

  void update(double t, double dt,
              const array::Array3D &u,
              const array::Array3D &v,
              const array::Scalar &ice_thickness,
              const array::Scalar &climatic_mass_balance,
              const array::Scalar &basal_melt_rate);

  const array::Array3D& layer_depths() const;

private:
  MaxTimestep max_timestep_impl(double t) const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  DiagnosticList diagnostics_impl() const;

  //! isochronal layer thicknesses
  std::shared_ptr<array::Array3D> m_layer_thickness;

  //! temporary storage needed for time stepping
  std::shared_ptr<array::Array3D> m_tmp;

  //! The index of the topmost isochronal layer.
  int m_top_layer;

  size_t m_time_index;
  std::vector<double> m_deposition_times;
};

} // end of namespace pism
