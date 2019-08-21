/* Copyright (C) 2016, 2017, 2019 PISM Authors
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
  GeometryEvolution(IceGrid::ConstPtr grid);
  ~GeometryEvolution();

  void init(const InputOptions &opts);

  void flow_step(const Geometry &ice_geometry, double dt,
                 const IceModelVec2V    &advective_velocity,
                 const IceModelVec2Stag &diffusive_flux,
                 const IceModelVec2Int  &velocity_bc_mask,
                 const IceModelVec2Int  &thickness_bc_mask);

  void source_term_step(const Geometry &geometry, double dt,
                        const IceModelVec2Int &thickness_bc_mask,
                        const IceModelVec2S   &surface_mass_flux,
                        const IceModelVec2S   &basal_melt_rate);

  void apply_flux_divergence(Geometry &geometry) const;
  void apply_mass_fluxes(Geometry &geometry) const;

  const IceModelVec2S& top_surface_mass_balance() const;
  const IceModelVec2S& bottom_surface_mass_balance() const;

  const IceModelVec2S& thickness_change_due_to_flow() const;
  const IceModelVec2S& area_specific_volume_change_due_to_flow() const;

  const IceModelVec2S& conservation_error() const;

  // diagnostic
  const IceModelVec2Stag& flux_staggered() const;
  const IceModelVec2S& flux_divergence() const;

  // "regional" setup
  virtual void set_no_model_mask(const IceModelVec2Int &mask);
protected:
  std::map<std::string,Diagnostic::Ptr> diagnostics_impl() const;

  virtual void init_impl(const InputOptions &opts);

  void update_in_place(double dt,
                       const IceModelVec2S& bed_elevation,
                       const IceModelVec2S& sea_level,
                       const IceModelVec2S& flux_divergence,
                       IceModelVec2S& ice_thickness,
                       IceModelVec2S& area_specific_volume);

  void residual_redistribution_iteration(const IceModelVec2S& bed_topography,
                                         const IceModelVec2S& sea_level,
                                         IceModelVec2S& ice_surface_elevation,
                                         IceModelVec2S& ice_thickness,
                                         IceModelVec2CellType& cell_type,
                                         IceModelVec2S& Href,
                                         IceModelVec2S& H_residual,
                                         bool &done);

  virtual void compute_interface_fluxes(const IceModelVec2CellType &cell_type,
                                        const IceModelVec2S        &ice_thickness,
                                        const IceModelVec2V        &velocity,
                                        const IceModelVec2Int      &velocity_bc_mask,
                                        const IceModelVec2Stag     &diffusive_flux,
                                        IceModelVec2Stag           &output);

  virtual void compute_flux_divergence(const IceModelVec2Stag &flux_staggered,
                                       const IceModelVec2Int &thickness_bc_mask,
                                       IceModelVec2S &flux_fivergence);

  virtual void ensure_nonnegativity(const IceModelVec2S &ice_thickness,
                                    const IceModelVec2S &area_specific_volume,
                                    IceModelVec2S &thickness_change,
                                    IceModelVec2S &area_specific_volume_change,
                                    IceModelVec2S &conservation_error);

  virtual void set_no_model_mask_impl(const IceModelVec2Int &mask);

  // note: cells with area_specific_volume > 0 do not experience changes due to surface and basal
  // mass balance sources
  virtual void compute_surface_and_basal_mass_balance(double dt,
                                                      const IceModelVec2Int      &thickness_bc_mask,
                                                      const IceModelVec2S        &ice_thickness,
                                                      const IceModelVec2CellType &cell_type,
                                                      const IceModelVec2S        &surface_mass_flux,
                                                      const IceModelVec2S        &basal_melt_rate,
                                                      IceModelVec2S              &effective_SMB,
                                                      IceModelVec2S              &effective_BMB);
protected:
  struct Impl;
  Impl *m_impl;
};

class RegionalGeometryEvolution : public GeometryEvolution {
public:
  RegionalGeometryEvolution(IceGrid::ConstPtr grid);

protected:
  void set_no_model_mask_impl(const IceModelVec2Int &mask);

  void compute_interface_fluxes(const IceModelVec2CellType &cell_type,
                                const IceModelVec2S        &ice_thickness,
                                const IceModelVec2V        &velocity,
                                const IceModelVec2Int      &velocity_bc_mask,
                                const IceModelVec2Stag     &diffusive_flux,
                                IceModelVec2Stag           &output);

  void compute_surface_and_basal_mass_balance(double dt,
                                              const IceModelVec2Int      &thickness_bc_mask,
                                              const IceModelVec2S        &ice_thickness,
                                              const IceModelVec2CellType &cell_type,
                                              const IceModelVec2S        &surface_mass_flux,
                                              const IceModelVec2S        &basal_melt_rate,
                                              IceModelVec2S              &effective_SMB,
                                              IceModelVec2S              &effective_BMB);
private:
  IceModelVec2Int m_no_model_mask;
};

} // end of namespace pism

#endif /* GEOMETRYEVOLUTION_H */
