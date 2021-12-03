/* Copyright (C) 2021 PISM Authors
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

#include "NoGLRetreat.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/pism_utilities.hh" // combine()

namespace pism {
namespace surface {

NoGLRetreat::NoGLRetreat(IceGrid::ConstPtr grid,
                         std::shared_ptr<SurfaceModel> input)
  : SurfaceModel(grid, input),
    m_smb_adjustment(grid, "smb_adjustment", WITHOUT_GHOSTS),
    m_min_ice_thickness(grid, "minimum_ice_thickness", WITHOUT_GHOSTS) {

  m_smb_adjustment.metadata()["units"] = "kg m-2 s-1";

  m_mass_flux    = allocate_mass_flux(grid);
  m_accumulation = allocate_accumulation(grid);
  m_melt         = allocate_melt(grid);
  m_runoff       = allocate_runoff(grid);
}

void NoGLRetreat::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing a SMB adjustment preventing grounding line retreat...\n");

  const auto &ice_thickness = geometry.ice_thickness;
  const auto &sea_level     = geometry.sea_level_elevation;
  const auto &bed           = geometry.bed_elevation;

  double rho_i = m_config->get_number("constants.ice.density");
  double rho_w = m_config->get_number("constants.sea_water.density");
  double eps = m_config->get_number("geometry.ice_free_thickness_standard");

  IceModelVec::AccessList list{&sea_level, &bed, &ice_thickness,
                               &m_min_ice_thickness};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double H_min = 0.0;
    if (sea_level(i, j) > bed(i, j)) {
      // compute the minimum grounded ice thickess
      H_min = (sea_level(i, j) - bed(i, j)) * (rho_w / rho_i) + eps;

      // exclude areas where the ice thickness is below this threshold (i.e. the area is
      // ice free or covered by floating ice)
      if (ice_thickness(i, j) < H_min) {
        H_min = 0.0;
      }
    }

    m_min_ice_thickness(i, j) = H_min;

  }
}

void NoGLRetreat::update_impl(const Geometry &geometry, double t, double dt) {

  m_input_model->update(geometry, t, dt);

  const auto &mass_flux     = m_input_model->mass_flux();
  const auto &ice_thickness = geometry.ice_thickness;

  double rho_i = m_config->get_number("constants.ice.density");

  IceModelVec::AccessList list{&mass_flux,
                               &ice_thickness,
                               m_mass_flux.get(),
                               &m_smb_adjustment,
                               &m_min_ice_thickness};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double SMB_old = mass_flux(i, j);
    double SMB_new = SMB_old;

    // convert from kg/(m^2 s) to m:
    double H_min = m_min_ice_thickness(i, j);
    double H     = ice_thickness(i, j);
    double dH    = mass_flux(i, j) * (dt / rho_i);
    double H_new = H + dH;

    if (H_min > 0.0 and H_new < H_min) {
      SMB_new = (H_min - H) * (rho_i / dt);
    }

    (*m_mass_flux)(i, j) = SMB_new;
    m_smb_adjustment(i, j) = SMB_new - SMB_old;
  }

  dummy_accumulation(*m_mass_flux, *m_accumulation);
  dummy_melt(*m_mass_flux, *m_melt);
  dummy_runoff(*m_mass_flux, *m_runoff);
}

const IceModelVec2S& NoGLRetreat::mass_flux_impl() const {
  return *m_mass_flux;
}

const IceModelVec2S& NoGLRetreat::accumulation_impl() const {
  return *m_accumulation;
}

const IceModelVec2S& NoGLRetreat::melt_impl() const {
  return *m_melt;
}

const IceModelVec2S& NoGLRetreat::runoff_impl() const {
  return *m_runoff;
}

const IceModelVec2S& NoGLRetreat::smb_adjustment() const {
  return m_smb_adjustment;
}

namespace diagnostics {

class SMBAdjustment : public DiagAverageRate<NoGLRetreat>
{
public:
  SMBAdjustment(const NoGLRetreat *m)
    : DiagAverageRate<NoGLRetreat>(m,
                                   "no_gl_retreat_smb_adjustment",
                                   RATE)
  {

    m_vars = {{m_sys, "no_gl_retreat_smb_adjustment"}};
    m_accumulator.metadata()["units"] = "kg m-2";

    set_attrs("SMB adjustment needed to maintain grounded ice extent",
              "",               // no standard name
              "kg m-2 s-1",
              "kg m-2 year-1",
              0);
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = {to_internal(m_fill_value)};
  }

protected:
  const IceModelVec2S& model_input() {
    return model->smb_adjustment();
  }
};

} // end of namespace diagnostics

DiagnosticList NoGLRetreat::diagnostics_impl() const {
  return combine({{"no_gl_retreat_smb_adjustment",
          Diagnostic::Ptr(new diagnostics::SMBAdjustment(this))}},
    m_input_model->diagnostics());
}

} // end of namespace surface
} // end of namespace pism
