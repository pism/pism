// Copyright (C) 2011, 2012 Constantine Khroulev
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

#include "iceModel.hh"
#include "IceGrid.hh"

#if (PISM_USE_PROJ4==1)

#include <proj_api.h>

//! Computes the area of a triangle using vector cross product.
static PetscReal triangle_area(PetscReal *A, PetscReal *B, PetscReal *C) {
  PetscReal V1[3], V2[3];
  for (int j = 0; j < 3; ++j) {
    V1[j] = B[j] - A[j];
    V2[j] = C[j] - A[j];
  }

  return 0.5*sqrt(PetscSqr(V1[1]*V2[2] - V2[1]*V1[2]) +
		  PetscSqr(V1[0]*V2[2] - V2[0]*V1[2]) +
		  PetscSqr(V1[0]*V2[1] - V2[0]*V1[1]));
}

PetscErrorCode IceModel::compute_cell_areas() {
  PetscErrorCode ierr;
  projPJ pism, lonlat, geocent;

  if (config.get_flag("correct_cell_areas") == false ||
      mapping.has("proj4") == false) {

    ierr = cell_area.set(grid.dx * grid.dy); CHKERRQ(ierr);

    return 0;
  }

  string proj_string = mapping.get_string("proj4");

  lonlat = pj_init_plus("+proj=latlong +datum=WGS84 +ellps=WGS84");
  if (lonlat == NULL) {
    PetscPrintf(grid.com, "PISM ERROR: projection initialization failed "
                "('+proj=latlong +datum=WGS84 +ellps=WGS84').\n");
    PISMEnd();
  }

  geocent = pj_init_plus("+proj=geocent +datum=WGS84 +ellps=WGS84");
  if (geocent == NULL) {
    PetscPrintf(grid.com, "PISM ERROR: projection initialization failed "
                "('+proj=geocent +datum=WGS84 +ellps=WGS84').\n");
    PISMEnd();
  }

  pism = pj_init_plus(proj_string.c_str());
  if (pism == NULL) {
    PetscPrintf(grid.com, "PISM ERROR: proj.4 string '%s' is invalid. Exiting...\n",
                proj_string.c_str());
    PISMEnd();
  }

  ierr = verbPrintf(2,grid.com,
		    "* Computing cell areas, latitude and longitude\n"
                    "  using projection parameters...\n"); CHKERRQ(ierr);

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

  double dx2 = 0.5 * grid.dx, dy2 = 0.5 * grid.dy;

  ierr = vLatitude.begin_access(); CHKERRQ(ierr);
  ierr = vLongitude.begin_access(); CHKERRQ(ierr);
  ierr =  cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      double x = grid.x[i], y = grid.y[j], Z;

      // compute the cell area:
      double x_nw = x - dx2, y_nw = y + dy2;
      Z = 0;
      pj_transform(pism, geocent, 1, 1, &x_nw, &y_nw, &Z);
      PetscReal nw[3] = {x_nw, y_nw, Z};

      double x_ne = x + dx2, y_ne = y + dy2;
      Z = 0;
      pj_transform(pism, geocent, 1, 1, &x_ne, &y_ne, &Z);
      PetscReal ne[3] = {x_ne, y_ne, Z};

      double x_se = x + dx2, y_se = y - dy2;
      Z = 0;
      pj_transform(pism, geocent, 1, 1, &x_se, &y_se, &Z);
      PetscReal se[3] = {x_se, y_se, Z};

      double x_sw = x - dx2, y_sw = y - dy2;
      Z = 0;
      pj_transform(pism, geocent, 1, 1, &x_sw, &y_sw, &Z);
      PetscReal sw[3] = {x_sw, y_sw, Z};

      cell_area(i, j) = triangle_area(sw, se, ne) + triangle_area(ne, nw, sw);

      // compute lon,lat coordinates:
      pj_transform(pism, lonlat, 1, 1, &x, &y, NULL);
      // NB! proj.4 converts x,y pairs into lon,lat pairs in *radians*.
      vLongitude(i, j) = x * RAD_TO_DEG;
      vLatitude(i, j)  = y * RAD_TO_DEG;
    }
  }
  ierr =  cell_area.end_access(); CHKERRQ(ierr);  
  ierr = vLongitude.end_access(); CHKERRQ(ierr);
  ierr = vLatitude.end_access(); CHKERRQ(ierr);

  pj_free(pism);
  pj_free(lonlat);
  pj_free(geocent);

  return 0;
}

#elif (PISM_USE_PROJ4==0)      // no proj.4

PetscErrorCode IceModel::compute_cell_areas() {
  PetscErrorCode ierr;

  // proj.4 was not found; use uncorrected areas.
  ierr = cell_area.set(grid.dx * grid.dy); CHKERRQ(ierr);

  return 0;
}

#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif
