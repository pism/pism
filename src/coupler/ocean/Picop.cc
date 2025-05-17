// Copyright (C) 2012-2019, 2021, 2022, 2023, 2024, 2025 Constantine Khrulev, Ricarda Winkelmann, Ronja Reese, Torsten
// Albrecht, and Matthias Mengel
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
//
// Please cite this model as:
// 1.
// Antarctic sub-shelf melt rates via PICO
// R. Reese, T. Albrecht, M. Mengel, X. Asay-Davis and R. Winkelmann
// The Cryosphere, 12, 1969-1985, (2018)
// DOI: 10.5194/tc-12-1969-2018
//
// 2.
// A box model of circulation and melting in ice shelf caverns
// D. Olbers & H. Hellmer
// Ocean Dynamics (2010), Volume 60, Issue 1, pp 141–153
// DOI: 10.1007/s10236-009-0252-z
//
// 3.
// PICOP, a new ocean melt parameterization under ice shelves
// combining PICO and a plume model.
// T. Pelle, M. Morlighem, J.H. Bondzio
// The Cryosphere, 13, 1043-49, (2019)
// DOI: 10.5194/tc-13-1043-2019

#include <cmath>

#include <gsl/gsl_math.h> // GSL_NAN

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/Scalar.hh"

#include "pism/coupler/ocean/Picop.hh"
#include "pism/coupler/ocean/PicoGeometry.hh"
#include "pism/coupler/ocean/PicoPhysics.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/Logger.hh"

namespace pism {

namespace stressbalance {
class StressBalance;
}

namespace ocean {

Picop::Picop(std::shared_ptr<const Grid> grid, std::shared_ptr<stressbalance::StressBalance> stressbalance):
    CompleteOceanModel(grid, std::shared_ptr<OceanModel>()),
    m_stress_balance(stressbalance),
    m_pico(std::make_shared<Pico>(grid)),
    m_grounding_line_elevation(grid, "picop_grounding_line_elevation"),
    m_theta_ocean(m_pico->get_temperature()),
    m_salinity_ocean(m_pico->get_salinity()),
    m_geometry(grid),
    m_velocity(grid, "ghosted_velocity")

{

  ForcingOptions opt(*m_grid->ctx(), "ocean.picop");

  m_grounding_line_elevation.metadata(0).long_name("grounding line elevation").units("m");
  m_grounding_line_elevation.metadata()["_FillValue"] = { 0.0 };
  m_grounding_line_elevation.set(0.0);

}

void Picop::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2, "* Initializing the Potsdam Ice-shelf Cavity mOdel / Plume for the ocean ...\n");

    
  compute_grounding_line_elevation(geometry, m_grounding_line_elevation);
}

void Picop::define_model_state_impl(const File &output) const {


  OceanModel::define_model_state_impl(output);
}

void Picop::write_model_state_impl(const File &output) const {


  OceanModel::write_model_state_impl(output);
}


void Picop::update_impl(const Geometry &geometry, double t, double dt) {

  (void) t;
  (void) dt;
  compute_grounding_line_elevation(geometry, m_grounding_line_elevation);
}


MaxTimestep Picop::max_timestep_impl(double t) const {
  (void) t;

  return MaxTimestep("ocean picop");
}


void Picop::compute_grounding_line_elevation(const Geometry &geometry,
                                             array::Scalar &grounding_line_elevation) const {
  const array::Scalar &bed        = geometry.bed_elevation;
  const array::Scalar &H          = geometry.ice_thickness;
  const array::Scalar &cell_type  = geometry.cell_type;
  const array::Scalar &z_s        = geometry.sea_level_elevation;

  const int Mx = m_grid->Mx();
  const int My = m_grid->My();
  const double dx = m_grid->dx();
  const double dy = m_grid->dy();

  array::AccessScope scope{&bed, &H, &z_s, &cell_type, &m_velocity, &grounding_line_elevation};

  // Step 1: Initialize zgl0 at grounding line: bed elevation
  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();
    grounding_line_elevation(i, j) = (cell_type(i, j) == MASK_GROUNDED) ? bed(i, j) : 0.0;
  }

  // Step 2: Advection-diffusion loop
  const double eps = 1e-8;
  const int max_iter = 100;
  const double tol = 1e-4;
  double residual = 0.0;

  using std::pow;
  using std::sqrt;

  for (int iter = 0; iter < max_iter; ++iter) {
    residual = 0.0;

    for (int j = 0; j < My; ++j) {
      for (int i = 0; i < Mx; ++i) {

        const double mag = sqrt(pow(m_velocity(i, j).u, 2.0) + pow(m_velocity(i, j).v, 2.0)) ;

        double u = m_velocity(i, j).u;
        double v = m_velocity(i, j).v;
        u /= mag;
        v /= mag;

        // x-derivative
        double dzdx;
        if (i > 0 && i < Mx - 1) {
          dzdx = (grounding_line_elevation(i + 1, j) - grounding_line_elevation(i - 1, j)) / (2.0 * dx);
        } else if (i == 0 && i + 2 < Mx) {
          dzdx = (-3 * grounding_line_elevation(i, j) + 4 * grounding_line_elevation(i + 1, j)
                  - grounding_line_elevation(i + 2, j)) / (2.0 * dx);
        } else if (i == Mx - 1 && i >= 2) {
          dzdx = (3 * grounding_line_elevation(i, j) - 4 * grounding_line_elevation(i - 1, j)
                  + grounding_line_elevation(i - 2, j)) / (2.0 * dx);
        } else {
          dzdx = 0.0; // fallback
        }

        // y-derivative
        double dzdy;
        if (j > 0 && j < My - 1) {
          dzdy = (grounding_line_elevation(i, j + 1) - grounding_line_elevation(i, j - 1)) / (2.0 * dy);
        } else if (j == 0 && j + 2 < My) {
          dzdy = (-3 * grounding_line_elevation(i, j) + 4 * grounding_line_elevation(i, j + 1)
                  - grounding_line_elevation(i, j + 2)) / (2.0 * dy);
        } else if (j == My - 1 && j >= 2) {
          dzdy = (3 * grounding_line_elevation(i, j) - 4 * grounding_line_elevation(i, j - 1)
                  + grounding_line_elevation(i, j - 2)) / (2.0 * dy);
        } else {
          dzdy = 0.0; // fallback
        }

        // Laplacian for diffusion
        double laplacian = 0.0;
        if (i > 0 && i < Mx - 1 && j > 0 && j < My - 1) {
          laplacian = (grounding_line_elevation(i + 1, j) + grounding_line_elevation(i - 1, j)
                     + grounding_line_elevation(i, j + 1) + grounding_line_elevation(i, j - 1)
                     - 4.0 * grounding_line_elevation(i, j)) / (dx * dx);
        }

        // Update
        double update = -(u * dzdx + v * dzdy + eps * laplacian);

        // Boundary condition: grounding_line_elevation == bed outside shelf
        if (cell_type(i, j) != MASK_FLOATING) {
          grounding_line_elevation(i, j) = bed(i,j);
        };
        grounding_line_elevation(i, j) += 0.1 * update;
        residual += update * update;
      }
    }

    residual = std::sqrt(GlobalSum(m_grid->com, residual));
    if (residual < tol) break;
  }

  // Step 3: Clip to base elevation
  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      double zb = z_s(i, j) - H(i, j);  // base of ice shelf
      if (grounding_line_elevation(i, j) > zb) {
        grounding_line_elevation(i, j) = zb;
      }
    }
  }

  m_log->message(2, "Computed grounding line elevation field with residual = %.3e\n", residual);
}

// Write diagnostic variables to extra files if requested
DiagnosticList Picop::diagnostics_impl() const {

  DiagnosticList result = {
    { "picop_grounding_line_elevation", Diagnostic::wrap(m_grounding_line_elevation) },
  };

  return combine(result, OceanModel::diagnostics_impl());
}


} // end of namespace ocean
} // end of namespace pism
