/* Copyright (C) 2025 PISM Authors
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

#ifndef PISM_GRIDINFO_H
#define PISM_GRIDINFO_H

#include <string>
#include <vector>

namespace pism {
namespace grid {

typedef enum {NOT_PERIODIC = 0, X_PERIODIC = 1, Y_PERIODIC = 2, XY_PERIODIC = 3} Periodicity;

typedef enum {CELL_CENTER, CELL_CORNER} Registration;

Registration string_to_registration(const std::string &keyword);
std::string registration_to_string(Registration registration);

Periodicity string_to_periodicity(const std::string &keyword);
std::string periodicity_to_string(Periodicity p);

class GridInfo {
public:
  //! x-coordinate of the domain center
  double x0;
  //! y-coordinate of the domain center
  double y0;
  //! domain half-width in the X direction
  double Lx;
  //! domain half-width in the Y direction
  double Ly;

  //! x coordinates
  std::vector<double> x;
  //! y coordinates
  std::vector<double> y;
};

class DistributedGridInfo : public GridInfo {
public:
  //! Grid periodicity in X and Y directions
  grid::Periodicity periodicity;
  //! Grid registration (cell center or cell corner)
  grid::Registration registration;

  //! Grid spacing in the X direction
  double dx;
  //! Grid spacing in the Y direction
  double dy;
  //! Cell area (meters^2) (same as dx*dy)
  double cell_area;

  //! Starting index (in the X direction) of the patch owned by the current MPI rank
  unsigned int xs;
  //! Number of grid points (in the X direction) of the patch owned by the current MPI
  //! rank
  unsigned int xm;
  //! Starting index of the patch owned by the current MPI rank (X direction)
  unsigned int ys;
  //! Number of grid points (in the Y direction) of the patch owned by the current MPI
  //! rank
  unsigned int ym;

  //! Total number of grid points in the X direction
  unsigned int Mx;
  //! Total number of grid points in the Y direction
  unsigned int My;

  //! Number of grid points in the largest patch (max(xm*ym) over all MPI ranks)
  int max_patch_size;

  //! Current MPI rank
  int rank;

  //! MPI Communicator size
  int size;
};

} // namespace grid
} // namespace pism

#endif /* PISM_GRIDINFO_H */
