// Copyright (C) 2012-2019, 2021, 2022, 2023, 2024, 2025, 2026 Constantine Khrulev, Ricarda Winkelmann, Ronja Reese, Torsten
// Albrecht, Matthias Mengel, and Andy Aschwanden
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef PISM_PICOP_H
#define PISM_PICOP_H

#include "pism/coupler/ocean/CompleteOceanModel.hh"
#include "pism/coupler/ocean/Pico.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {

namespace ocean {

class PicopPhysics;

//! Buoyant-plume sub-shelf melt (Lazeroms et al. 2018, Pelle et al. 2019) with the
//! subglacial-discharge extension of Pelle et al. (2023), driven by ambient ocean
//! conditions supplied by a derived class.
/*!
 * Holds everything the plume needs that does not depend on where the ambient
 * temperature and salinity come from: the grounding-line elevation transported along
 * the flow, the shelf-base elevation and slope, the subglacial-discharge field and the
 * melt-rate computation itself. Picop feeds it the PICO box-model output; Plume feeds it
 * fields read from a file.
 */
class PlumeModel : public CompleteOceanModel {
public:
  virtual ~PlumeModel() = default;

protected:
  PlumeModel(std::shared_ptr<const Grid> g);

  //! Update the shelf-base elevation, the grounding-line elevation and the local slope.
  void update_geometry(const Inputs &inputs);

  //! Compute the plume melt rate from ambient temperature `T_a` (kelvin) and salinity
  //! `S_a` (g/kg) on floating cells and set the shelf-base mass flux and the average
  //! water-column pressure. Call after update_geometry().
  void update_melt_rate(const Inputs &inputs, const PicopPhysics &physics,
                        const array::Scalar &T_a, const array::Scalar &S_a);

  //! PICOP's spatial diagnostics (`picop_*`), with `T_a` and `S_a` as
  //! `picop_temperature` and `picop_salinity`.
  DiagnosticList plume_diagnostics(const array::Scalar &T_a, const array::Scalar &S_a) const;

  array::Scalar1 m_basal_melt_rate;
  array::Scalar1 m_grounding_line_elevation;
  array::Scalar1 m_shelf_base_elevation;
  array::Scalar1 m_local_slope;
  array::Scalar1 m_fresh_water_melt_rate;
  //! subglacial discharge flux q_sg(x,y) on floating cells (m^2 s^-1)
  array::Scalar m_discharge_flux;

  //! transported tracers for along-flow discharge distribution: source discharge q_sg0,
  //! source governing length scale 5L', and along-flow path distance from the source
  array::Scalar1 m_disch_q0, m_disch_L5, m_disch_s;

  //! whether to add the fresh-water (subglacial discharge) melt contribution
  bool m_add_fresh_water_melt;

  //! how the subglacial discharge plume q_sg(x,y) is distributed onto floating cells
  enum DischargeMethod { DISCHARGE_ISOTROPIC, DISCHARGE_DOWNSTREAM_GATE, DISCHARGE_ALONG_FLOW };
  DischargeMethod m_discharge_method;

  array::Vector m_flow_direction;
  array::Scalar m_work;

  //! temporary storage for the shelf base gradient
  array::Staggered1 m_zb_x, m_zb_y;

private:
  void compute_melt_rate(const Inputs &inputs,
                         const PicopPhysics &physics,
                         const array::Scalar &T_a,
                         const array::Scalar &S_a,
                         array::Scalar1 &result);

  void compute_grounding_line_elevation(const Inputs &inputs,
                                        array::Scalar1 &result);

  void compute_shelf_base_elevation(const Inputs &inputs,
                                        array::Scalar1 &result);

  void compute_local_slope(const Inputs &inputs,
                                        array::Scalar1 &result);

  //! Build q_sg(x,y) on floating cells from grounding-line discharge outflows.
  void build_discharge_field(const Inputs &inputs,
                             const PicopPhysics &physics,
                             const array::Scalar &T_a,
                             const array::Scalar &S_a);
};

//! Implements the PICO ocean model published in The Cryosphere (2018)
//! and adds the Plume model Pelle et al (2019).
//!
class Picop : public PlumeModel {
public:
  Picop(std::shared_ptr<const Grid> g);
  virtual ~Picop() = default;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Inputs &inputs, double t, double dt);
  MaxTimestep max_timestep_impl(double t, const CFLData *cfl_data) const;

  std::set<VariableMetadata> state_impl() const;
  void write_state_impl(const OutputFile &output) const;

  std::map<std::string, Diagnostic::Ptr> spatial_diagnostics_impl() const;

private:
  std::shared_ptr<Pico> m_pico;

  //! ambient temperature and salinity: PICO's box values on floating cells
  const array::Scalar &m_theta_ocean;
  const array::Scalar &m_salinity_ocean;
};
} // end of namespace ocean
} // end of namespace pism

#endif /* PISM_PICOP_H */
