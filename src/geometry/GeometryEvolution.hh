/* Copyright (C) 2016, 2017, 2019, 2020, 2022 PISM Authors
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

#ifndef GEOMETRYEVOLUTION_H
#define GEOMETRYEVOLUTION_H

#include "Geometry.hh"
#include "pism/util/Component.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {

/*!
 * NB! Write this in a way that does not use ghosts of input fields (copy to temp. storage and
 * communicate).
 *
 * The promise:
 *
 * H_change + Href_change = dt * (SMB_rate + BMB_rate - flux_divergence) + conservation_error
 *
 * Href == 0 if H > 0
 */
class GeometryEvolution : public Component {
public:
  GeometryEvolution(std::shared_ptr<const Grid> grid);
  ~GeometryEvolution();

  void init(const InputOptions &opts);

  void reset();

  void flow_step(const Geometry &ice_geometry, double dt,
                 const array::Vector    &advective_velocity,
                 const array::Staggered &diffusive_flux,
                 const array::Scalar  &thickness_bc_mask);

  void source_term_step(const Geometry &geometry, double dt,
                        const array::Scalar &thickness_bc_mask,
                        const array::Scalar   &surface_mass_flux,
                        const array::Scalar   &basal_melt_rate);

  void apply_flux_divergence(Geometry &geometry) const;
  void apply_mass_fluxes(Geometry &geometry) const;

  const array::Scalar& top_surface_mass_balance() const;
  const array::Scalar& bottom_surface_mass_balance() const;

  const array::Scalar& thickness_change_due_to_flow() const;
  const array::Scalar& area_specific_volume_change_due_to_flow() const;

  const array::Scalar& conservation_error() const;

  // diagnostic
  const array::Staggered1& flux_staggered() const;
  const array::Scalar& flux_divergence() const;

  // "regional" setup
  virtual void set_no_model_mask(const array::Scalar &mask);
protected:
  std::map<std::string,Diagnostic::Ptr> diagnostics_impl() const;

  virtual void init_impl(const InputOptions &opts);

  void update_in_place(double dt,
                       const array::Scalar& bed_elevation,
                       const array::Scalar& sea_level,
                       const array::Scalar& flux_divergence,
                       array::Scalar& ice_thickness,
                       array::Scalar& area_specific_volume);

  void residual_redistribution_iteration(const array::Scalar &bed_topography,
                                         const array::Scalar &sea_level,
                                         array::Scalar1      &ice_surface_elevation,
                                         array::Scalar       &ice_thickness,
                                         array::CellType1    &cell_type,
                                         array::Scalar       &Href,
                                         array::Scalar       &H_residual,
                                         bool                &done);

  virtual void compute_interface_fluxes(const array::CellType1 &cell_type,
                                        const array::Scalar        &ice_thickness,
                                        const array::Vector        &velocity,
                                        const array::Staggered     &diffusive_flux,
                                        array::Staggered           &output);

  virtual void compute_flux_divergence(double dt,
                                       const array::Staggered1 &flux_staggered,
                                       const array::Scalar &thickness_bc_mask,
                                       array::Scalar &conservation_error,
                                       array::Scalar &flux_fivergence);

  virtual void ensure_nonnegativity(const array::Scalar &ice_thickness,
                                    const array::Scalar &area_specific_volume,
                                    array::Scalar &thickness_change,
                                    array::Scalar &area_specific_volume_change,
                                    array::Scalar &conservation_error);

  virtual void set_no_model_mask_impl(const array::Scalar &mask);

  // note: cells with area_specific_volume > 0 do not experience changes due to surface and basal
  // mass balance sources
  virtual void compute_surface_and_basal_mass_balance(double dt,
                                                      const array::Scalar      &thickness_bc_mask,
                                                      const array::Scalar        &ice_thickness,
                                                      const array::CellType &cell_type,
                                                      const array::Scalar        &surface_mass_flux,
                                                      const array::Scalar        &basal_melt_rate,
                                                      array::Scalar              &effective_SMB,
                                                      array::Scalar              &effective_BMB);
protected:
  struct Impl;
  Impl *m_impl;
};

class RegionalGeometryEvolution : public GeometryEvolution {
public:
  RegionalGeometryEvolution(std::shared_ptr<const Grid> grid);

protected:
  void set_no_model_mask_impl(const array::Scalar &mask);

  void compute_interface_fluxes(const array::CellType1 &cell_type,
                                const array::Scalar        &ice_thickness,
                                const array::Vector        &velocity,
                                const array::Staggered     &diffusive_flux,
                                array::Staggered           &output);

  void compute_surface_and_basal_mass_balance(double dt,
                                              const array::Scalar      &thickness_bc_mask,
                                              const array::Scalar        &ice_thickness,
                                              const array::CellType &cell_type,
                                              const array::Scalar        &surface_mass_flux,
                                              const array::Scalar        &basal_melt_rate,
                                              array::Scalar              &effective_SMB,
                                              array::Scalar              &effective_BMB);
private:
  array::Scalar1 m_no_model_mask;
};

/*!
 * Compute the grounding line flux.
 *
 * The units of `result` are "kg m-2". Negative flux corresponds to ice moving into
 * the ocean, i.e. from grounded to floating areas.
 *
 * This convention makes it easier to compare this quantity to the surface mass balance or
 * calving fluxes.
 */
void grounding_line_flux(const array::CellType1 &cell_type,
                         const array::Staggered1 &flux,
                         double dt,
                         bool add_values,
                         array::Scalar &result);

double total_grounding_line_flux(const array::CellType1 &cell_type,
                                 const array::Staggered1 &flux,
                                 double dt);
} // end of namespace pism

#endif /* GEOMETRYEVOLUTION_H */
