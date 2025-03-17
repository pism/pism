/* Copyright (C) 2020, 2021, 2023, 2024, 2025 PISM Authors
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

#include "pism/stressbalance/blatter/ismip-hom/BlatterISMIPHOM.hh"

#include "pism/util/node_types.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace stressbalance {

static double A_surface(double x, double y, double L) {
  (void) y;
  (void) L;
  double alpha = 0.5 * (M_PI / 180.0); // 0.5 degrees
  return -x * tan(alpha);
}

static double A_bed(double x, double y, double L) {
  double omega = 2 * M_PI / L;
  return A_surface(x, y, L) - 1000.0 + 500.0 * sin(omega * x) * sin(omega * y);
}

static double B_bed(double x, double y, double L) {
  (void) y;
  double omega = 2 * M_PI / L;
  return A_surface(x, y, L) - 1000.0 + 500.0 * sin(omega * x);
}

static double C_surface(double x, double y, double L) {
  (void) y;
  (void) L;
  double alpha = 0.1 * (M_PI / 180.0); // 0.1 degrees
  return -x * tan(alpha);
}

static double C_bed(double x, double y, double L) {
  return C_surface(x, y, L) - 1000.0;
}

BlatterISMIPHOM::BlatterISMIPHOM(std::shared_ptr<const Grid> grid,
                                 int Mz,
                                 int coarsening_factor,
                                 ISMIPHOMTest test)
  : Blatter(grid, Mz, coarsening_factor),
    m_test(test),
    m_L(2.0 * grid->Lx()) {

  char testname[] = "ABCD";
  m_log->message(2, "Running ISMIP-HOM Experiment %c (L = %d km)...\n",
                 testname[m_test],
                 (int)(m_L * 1e-3));

  switch (m_test) {
  case HOM_A:
    m_b = A_bed;
    m_s = A_surface;
    break;
  case HOM_B:
    m_b = B_bed;
    // test B surface is the same as for test A
    m_s = A_surface;
    break;
  case HOM_C:
  case HOM_D:
    // test D geometry is the same as for test C
    m_b = C_bed;
    m_s = C_surface;
    break;
  default:
    m_b = nullptr;
    m_s = nullptr;
    break;
  }
}

void BlatterISMIPHOM::nodal_parameter_values(const fem::Q1Element3 &element,
                                             Parameters **P,
                                             int i,
                                             int j,
                                             int *node_type,
                                             double *bottom_elevation,
                                             double *ice_thickness,
                                             double *surface_elevation,
                                             double *sea_level) const {
  (void) P;

  // This method is called before we get a chance to "reset" the current element, so the
  // element argument does not yet know about its physical coordinates. This is why we
  // have to compute x and y "by hand".
  double
    x_min = m_grid->x0() - m_grid->Lx(),
    y_min = m_grid->y0() - m_grid->Ly(),
    dx    = m_grid->dx(),
    dy    = m_grid->dy();

  for (int n = 0; n < fem::q13d::n_chi; ++n) {
    auto I = element.local_to_global(i, j, 0, n);

    double
      x = x_min + I.i * dx,
      y = y_min + I.j * dy;

    node_type[n]        = NODE_INTERIOR;
    bottom_elevation[n] = m_b(x, y, m_L);
    ice_thickness[n]    = m_s(x, y, m_L) - m_b(x, y, m_L);

    if (surface_elevation != nullptr) {
      surface_elevation[n] = m_s(x, y, m_L);
    }

    if (sea_level != nullptr) {
      // set sea level low enough to ensure that all ice is grounded
      sea_level[n] = bottom_elevation[n] - 1.0;
    }
  }
}


} // end of namespace stressbalance
} // end of namespace pism
