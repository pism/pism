// Copyright (C) 2004-2024 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef PISM_GRID_H
#define PISM_GRID_H

#include "io/File.hh"
#include <cassert>
#include <memory> // shared_ptr
#include <string>
#include <vector>
#include <map>

#include <mpi.h>                // MPI_Comm

#include "pism/util/Interpolation1D.hh"

namespace pism {

class Config;
class Context;
class File;
class InputInterpolation;
class Logger;
class MappingInfo;
class Vars;

namespace petsc {
class DM;
} // end of namespace petsc

namespace units {
class System;
}

namespace grid {

typedef enum {UNKNOWN = 0, EQUAL, QUADRATIC} VerticalSpacing;
typedef enum {NOT_PERIODIC = 0, X_PERIODIC = 1, Y_PERIODIC = 2, XY_PERIODIC = 3} Periodicity;

typedef enum {CELL_CENTER, CELL_CORNER} Registration;

Registration string_to_registration(const std::string &keyword);
std::string registration_to_string(Registration registration);

Periodicity string_to_periodicity(const std::string &keyword);
std::string periodicity_to_string(Periodicity p);

VerticalSpacing string_to_spacing(const std::string &keyword);
std::string spacing_to_string(VerticalSpacing s);

//! @brief Contains parameters of an input file grid.
class InputGridInfo {
public:
  InputGridInfo(const File &file, const std::string &variable,
                std::shared_ptr<units::System> unit_system, Registration registration);

  void report(const Logger &log, int threshold, std::shared_ptr<units::System> s) const;

  // dimension lengths
  unsigned int t_len;
  //! x-coordinate of the domain center
  double x0;
  //! y-coordinate of the domain center
  double y0;
  //! domain half-width
  double Lx;
  //! domain half-height
  double Ly;
  //! minimal value of the z dimension
  double z_min;
  //! maximal value of the z dimension
  double z_max;

  //! x coordinates
  std::vector<double> x;
  //! y coordinates
  std::vector<double> y;
  //! z coordinates
  std::vector<double> z;

  std::string filename;

  /*!
   * Variable name *found in the file*, which may not match the one expected by PISM.
   */
  std::string variable_name;

  std::map<std::string, AxisType> dimension_types;

  bool longitude_latitude;
private:
  void reset();
};

//! Grid parameters; used to collect defaults before an Grid is allocated.
/* Make sure that all of
   - `horizontal_size_from_options()`
   - `horizontal_extent_from_options()`
   - `vertical_grid_from_options()`
   - `ownership_ranges_from_options()`

   are called *or* all data members (`Lx`, `Ly`, `x0`, `y0`, `Mx`, `My`, `z`, `periodicity`,
   `procs_x`, `procs_y`) are set manually before using an instance of GridParameters.

   Call `validate()` to check current parameters.
*/
class Parameters {
public:
  //! Initialize grid defaults from a configuration database.
  Parameters(const Config &config);

  //! Initialize grid defaults from a configuration database, but set Mx and My explicitly.
  Parameters(const Config &config, unsigned Mx_, unsigned My_, double Lx, double Ly);

  //! Initialize grid defaults from a NetCDF variable.
  Parameters(std::shared_ptr<units::System> unit_system, const File &file,
             const std::string &variable_name, Registration r);

  static Parameters FromGridDefinition(std::shared_ptr<units::System> unit_system, const File &file,
                                       const std::string &variable_name, Registration registration);

  //! Process -Lx, -Ly, -x0, -y0; set Lx, Ly, x0, y0.
  void horizontal_size_and_extent_from_options(const Config &config);
  //! Process -Mz and -Lz; set z;
  void vertical_grid_from_options(const Config &config);
  //! Re-compute ownership ranges. Uses current values of Mx and My.
  void ownership_ranges_from_options(const Config &config, unsigned int size);

  //! Validate data members.
  void validate() const;

  //! Domain half-width in the X direction
  double Lx;
  //! Domain half-width in the Y direction
  double Ly;
  //! Domain center in the X direction
  double x0;
  //! Domain center in the Y direction
  double y0;
  //! Number of grid points in the X direction
  unsigned int Mx;
  //! Number of grid points in the Y direction
  unsigned int My;
  //! Grid registration
  Registration registration;
  //! Grid periodicity
  Periodicity periodicity;
  //! Vertical levels
  std::vector<double> z;
  //! Processor ownership ranges in the X direction
  std::vector<unsigned int> procs_x;
  //! Processor ownership ranges in the Y direction
  std::vector<unsigned int> procs_y;

  //! Name of the variable used to initialize the instance (empty if not used)
  std::string variable_name;
private:
  Parameters() = default;
};
} // namespace grid

class Grid;

/** Iterator class for traversing the grid, including ghost points.
 *
 * Usage:
 *
 * `for (PointsWithGhosts p(grid, stencil_width); p; p.next()) { ... }`
 */
class PointsWithGhosts {
public:
  PointsWithGhosts(const Grid &grid, unsigned int stencil_width = 1);

  int i() const {
    return m_i;
  }
  int j() const {
    return m_j;
  }

  void next() {
    assert(not m_done);
    m_i += 1;
    if (m_i > m_i_last) {
      m_i = m_i_first;        // wrap around
      m_j += 1;
    }
    if (m_j > m_j_last) {
      m_j = m_j_first;        // ensure that indexes are valid
      m_done = true;
    }
  }

  operator bool() const {
    return not m_done;
  }
private:
  int m_i, m_j;
  int m_i_first, m_i_last, m_j_first, m_j_last;
  bool m_done;
};

/** Iterator class for traversing the grid (without ghost points).
 *
 * Usage:
 *
 * `for (Points p(grid); p; p.next()) { int i = p.i(), j = p.j(); ... }`
 */
class Points : public PointsWithGhosts {
public:
  Points(const Grid &g) : PointsWithGhosts(g, 0) {}
};


//! Describes the PISM grid and the distribution of data across processors.
/*!
  This class holds parameters describing the grid, including the vertical
  spacing and which part of the horizontal grid is owned by the processor. It
  contains the dimensions of the PISM (4-dimensional, x*y*z*time) computational
  box. The vertical spacing can be quite arbitrary.

  It creates and destroys a two dimensional `PETSc` `DA` (distributed array).
  The creation of this `DA` is the point at which PISM gets distributed across
  multiple processors.

  \section computational_grid Organization of PISM's computational grid

  PISM uses the class Grid to manage computational grids and their
  parameters.

  Computational grids PISM can use are
  - rectangular,
  - equally spaced in the horizintal (X and Y) directions,
  - distributed across processors in horizontal dimensions only (every column
  is stored on one processor only),
  - are periodic in both X and Y directions (in the topological sence).

  Each processor "owns" a rectangular patch of `xm` times `ym` grid points with
  indices starting at `xs` and `ys` in the X and Y directions respectively.

  The typical code performing a point-wise computation will look like

  \code
  for (Points p(grid); p; p.next()) {
  const int i = p.i(), j = p.j();
  field(i,j) = value;
  }
  \endcode

  For finite difference (and some other) computations we often need to know
  values at map-plane neighbors of a grid point. 

  We say that a patch owned by a processor is surrounded by a strip of "ghost"
  grid points belonging to patches next to the one in question. This lets us to
  access (read) values at all the eight neighbors of a grid point for *all*
  the grid points, including ones at an edge of a processor patch *and* at an
  edge of a computational domain.

  All the values *written* to ghost points will be lost next time ghost values
  are updated.

  Sometimes it is beneficial to update ghost values locally (for instance when
  a computation A uses finite differences to compute derivatives of a quantity
  produced using a purely local (point-wise) computation B). In this case the
  loop above can be modified to look like

  \code
  for (PointsWithGhosts p(grid, ghost_width); p; p.next()) {
  const int i = p.i(), j = p.j();
  field(i,j) = value;
  }
  \endcode
*/
class Grid {
public:
  ~Grid();

  Grid(std::shared_ptr<const Context> context, const grid::Parameters &p);

  static std::shared_ptr<Grid> Shallow(std::shared_ptr<const Context> ctx, double Lx, double Ly,
                                       double x0, double y0, unsigned int Mx, unsigned int My,
                                       grid::Registration r, grid::Periodicity p);

  static std::shared_ptr<Grid> FromFile(std::shared_ptr<const Context> ctx,
                                        const File &file,
                                        const std::vector<std::string> &var_names,
                                        grid::Registration r);

  static std::shared_ptr<Grid> FromOptions(std::shared_ptr<const Context> ctx);

  std::shared_ptr<petsc::DM> get_dm(unsigned int dm_dof, unsigned int stencil_width) const;

  std::shared_ptr<InputInterpolation> get_interpolation(const std::vector<double> &levels,
                                                        const File &input_file,
                                                        const std::string &variable_name,
                                                        InterpolationType type) const;

  void forget_interpolations();

  void report_parameters() const;

  void compute_point_neighbors(double X, double Y,
                               int &i_left, int &i_right,
                               int &j_bottom, int &j_top) const;
  std::vector<int> point_neighbors(double X, double Y) const;
  std::vector<double> interpolation_weights(double x, double y) const;

  unsigned int kBelowHeight(double height) const;

  int max_patch_size() const;

  std::shared_ptr<const Context> ctx() const;

  int xs() const;
  int xm() const;
  int ys() const;
  int ym() const;

  const std::vector<double>& x() const;
  double x(size_t i) const;

  const std::vector<double>& y() const;
  double y(size_t i) const;

  const std::vector<double>& z() const;
  double z(size_t i) const;

  double dx() const;
  double dy() const;
  double cell_area() const;

  unsigned int Mx() const;
  unsigned int My() const;
  unsigned int Mz() const;

  double Lx() const;
  double Ly() const;
  double Lz() const;
  double x0() const;
  double y0() const;

  const MappingInfo& get_mapping_info() const;
  void set_mapping_info(const MappingInfo &info);

  double dz_min() const;
  double dz_max() const;

  grid::Periodicity periodicity() const;
  grid::Registration registration() const;

  unsigned int size() const;
  int rank() const;

  const MPI_Comm com;

  Vars& variables();
  const Vars& variables() const;

  int pio_io_decomposition(int dof, int output_datatype) const;

  PointsWithGhosts points(unsigned int stencil_width = 0) const {
    return {*this, stencil_width};
  }

private:
  struct Impl;
  Impl *m_impl;

  // Hide copy constructor / assignment operator.
  Grid(const Grid &);
  Grid & operator=(const Grid &);
};

namespace grid {

std::vector<double> compute_vertical_levels(double new_Lz, unsigned int new_Mz,
                                            grid::VerticalSpacing spacing, double Lambda = 0.0);

//! Returns the distance from the point (i,j) to the origin.
double radius(const Grid &grid, int i, int j);

//! @brief Check if a point `(i,j)` is in the strip of `stripwidth` meters around the edge
//! of the computational domain.
inline bool in_null_strip(const Grid& grid, int i, int j, double strip_width) {
  return (strip_width >  0.0                               &&
          (grid.x(i)  <= grid.x(0) + strip_width           ||
           grid.x(i)  >= grid.x(grid.Mx()-1) - strip_width ||
           grid.y(j)  <= grid.y(0) + strip_width           ||
           grid.y(j)  >= grid.y(grid.My()-1) - strip_width));
}

/*!
 * Return `true` if a point `(i, j)` is at an edge of the domain defined by the `grid`.
 */
inline bool domain_edge(const Grid &grid, int i, int j) {
  return ((j == 0) or (j == (int)grid.My() - 1) or (i == 0) or (i == (int)grid.Mx() - 1));
}

std::array<unsigned, 2> nprocs(unsigned int size, unsigned int Mx,
                                       unsigned int My);

std::vector<unsigned int> ownership_ranges(unsigned int Mx, unsigned int Nx);

} // namespace grid

} // end of namespace pism

#endif  /* PISM_GRID_H */
