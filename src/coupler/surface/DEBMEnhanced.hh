// Copyright (C) 2026 PISM Authors
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

#ifndef PISM_DEBM_ENHANCED_H
#define PISM_DEBM_ENHANCED_H

#include <memory>
#include <string>
#include <vector>

#include "pism/coupler/surface/DEBMSimple.hh"

namespace pism {

namespace array {
class Forcing;
class AccessScope;
class Scalar;
} // namespace array

namespace surface {

class TerrainInsolation;

//! @brief dEBM-enhanced: dEBM-simple with the insolation-driven melt computed from a
//! prescribed daily surface-insolation field rather than the analytic top-of-atmosphere
//! parameterization.
/*!
 * dEBM-enhanced reuses all of dEBM-simple (temperature- and offset-driven melt, albedo,
 * refreezing, snow bookkeeping, diagnostics). The only difference is the source of the
 * insolation energy that drives the insolation-melt term: instead of the analytic
 * top-of-atmosphere insolation, it reads a daily, terrain-shaded surface-insolation field
 * (variable `insolation`, units J m-2) precomputed by `util/pism_compute_insolation`.
 *
 * All physical parameters are shared with dEBM-simple (the `surface.debm_simple.*`
 * configuration namespace); only the input file is configured under
 * `surface.debm_enhanced.*`.
 */
class DEBMEnhanced : public DEBMSimple {
public:
  DEBMEnhanced(std::shared_ptr<const Grid> g,
               std::shared_ptr<atmosphere::AtmosphereModel> input);
  virtual ~DEBMEnhanced();

protected:
  void init_impl(const Geometry &geometry) override;

  void update_insolation_input(double t, double dt,
                               const std::vector<double> &ts,
                               array::AccessScope &list) override;

  void insolation_energy_series(int i, int j,
                                const std::vector<DEBMSimpleOrbitalParameters> &orbital,
                                const std::vector<double> &ts,
                                double dt_sub,
                                double latitude,
                                std::vector<double> &result) const override;

  DiagnosticList spatial_diagnostics_impl() const override;

private:
  //! true when no input file is given and the insolation is computed internally
  bool m_compute_internally;

  //! prescribed daily surface insolation energy (J m-2), read from a file (file path)
  std::shared_ptr<array::Forcing> m_input_insolation;

  std::string m_insolation_file;
  bool m_insolation_periodic;

  //! terrain-horizon + shaded-insolation engine (internal-compute path)
  std::unique_ptr<TerrainInsolation> m_terrain;

  //! daily surface insolation energy (J m-2) computed internally for the current update
  std::shared_ptr<array::Scalar> m_computed_insolation;

  //! per-cell latitude (degrees north), cached from the geometry at initialization
  const array::Scalar *m_latitude;

  //! ice surface elevation, cached from the geometry at initialization (the horizon is
  //! recomputed from this evolving field every m_update_interval)
  const array::Scalar *m_surface_elevation;

  //! interval (seconds) between recomputations of the terrain horizon
  double m_update_interval;

  //! model time (seconds) of the most recent horizon computation
  double m_t_last_horizon;

  //! seconds in a day, used to convert daily energy to per-sub-step energy
  double m_seconds_per_day;
};

} // end of namespace surface
} // end of namespace pism

#endif /* PISM_DEBM_ENHANCED_H */
