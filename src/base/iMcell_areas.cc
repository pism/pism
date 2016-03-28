// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
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

#if (PISM_USE_PROJ4==1)
#include <proj_api.h>
#endif

#include "iceModel.hh"

#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"

namespace pism {

#if (PISM_USE_PROJ4==1)

//! Computes the area of a triangle using vector cross product.
static double triangle_area(double *A, double *B, double *C) {
  double V1[3], V2[3];
  for (int j = 0; j < 3; ++j) {
    V1[j] = B[j] - A[j];
    V2[j] = C[j] - A[j];
  }

  return 0.5*sqrt(PetscSqr(V1[1]*V2[2] - V2[1]*V1[2]) +
                  PetscSqr(V1[0]*V2[2] - V2[0]*V1[2]) +
                  PetscSqr(V1[0]*V2[1] - V2[0]*V1[1]));
}

void IceModel::compute_cell_areas() {
  projPJ pism, lonlat, geocent;

  if (not m_config->get_boolean("correct_cell_areas") ||
      not m_output_global_attributes.has_attribute("proj4")) {

    m_log->message(2,
                   "* Computing cell areas using grid spacing (dx = %f m, dy = %f m)...\n",
                   m_grid->dx(), m_grid->dy());

    m_cell_area.set(m_grid->dx() * m_grid->dy());

    return;
  }

  std::string proj_string = m_output_global_attributes.get_string("proj4");

  lonlat = pj_init_plus("+proj=latlong +datum=WGS84 +ellps=WGS84");
  if (lonlat == NULL) {
    throw RuntimeError("projection initialization failed\n"
                       "('+proj=latlong +datum=WGS84 +ellps=WGS84').");
  }

  geocent = pj_init_plus("+proj=geocent +datum=WGS84 +ellps=WGS84");
  if (geocent == NULL) {
    throw RuntimeError("projection initialization failed\n"
                       "('+proj=geocent +datum=WGS84 +ellps=WGS84').");
  }

  pism = pj_init_plus(proj_string.c_str());
  if (pism == NULL) {
    throw RuntimeError::formatted("proj.4 string '%s' is invalid.",
                                  proj_string.c_str());
  }

  m_log->message(2,
                 "* Computing cell areas, latitude and longitude\n"
                 "  using projection parameters (%s)...\n",
                 proj_string.c_str());

// Cell layout:
// (nw)        (ne)
// +-----------+
// |           |
// |           |
// |     o     |
// |           |
// |           |
// +-----------+
// (sw)        (se)

  double dx2 = 0.5 * m_grid->dx(), dy2 = 0.5 * m_grid->dy();

  IceModelVec::AccessList list;
  list.add(vLatitude);
  list.add(vLongitude);
  list.add(m_cell_area);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double x = m_grid->x(i), y = m_grid->y(j), Z;

    // compute the cell area:
    double x_nw = x - dx2, y_nw = y + dy2;
    Z = 0;
    pj_transform(pism, geocent, 1, 1, &x_nw, &y_nw, &Z);
    double nw[3] = {x_nw, y_nw, Z};

    double x_ne = x + dx2, y_ne = y + dy2;
    Z = 0;
    pj_transform(pism, geocent, 1, 1, &x_ne, &y_ne, &Z);
    double ne[3] = {x_ne, y_ne, Z};

    double x_se = x + dx2, y_se = y - dy2;
    Z = 0;
    pj_transform(pism, geocent, 1, 1, &x_se, &y_se, &Z);
    double se[3] = {x_se, y_se, Z};

    double x_sw = x - dx2, y_sw = y - dy2;
    Z = 0;
    pj_transform(pism, geocent, 1, 1, &x_sw, &y_sw, &Z);
    double sw[3] = {x_sw, y_sw, Z};

    m_cell_area(i, j) = triangle_area(sw, se, ne) + triangle_area(ne, nw, sw);

    // compute lon,lat coordinates:
    pj_transform(pism, lonlat, 1, 1, &x, &y, NULL);
    // NB! proj.4 converts x,y pairs into lon,lat pairs in *radians*.
    vLongitude(i, j) = x * RAD_TO_DEG;
    vLatitude(i, j)  = y * RAD_TO_DEG;
  }

  pj_free(pism);
  pj_free(lonlat);
  pj_free(geocent);
}

#elif (PISM_USE_PROJ4==0)      // no proj.4

void IceModel::compute_cell_areas() {
  // proj.4 was not found; use uncorrected areas.
  m_cell_area.set(m_grid->dx() * m_grid->dy());
}

#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

} // end of namespace pism
