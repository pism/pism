// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __grid_hh
#define __grid_hh

#include <petscdmda.h>
#include <vector>
#include <string>
#include <map>
#include "PISMUnits.hh"

#ifdef PISM_USE_TR1
#include <tr1/memory>
#else
#include <memory>
#endif

#include <cassert>

namespace pism {

class Time;
class Prof;
class Config;

typedef enum {UNKNOWN = 0, EQUAL, QUADRATIC} SpacingType;
typedef enum {NONE = 0, NOT_PERIODIC =0, X_PERIODIC = 1, Y_PERIODIC = 2, XY_PERIODIC = 3} Periodicity;

/** Wrapper around PETSc's DM. Simplifies memory management.
 *
 * The constructor takes ownership of the dm argument passed to it.
 *
 * The destructor call DMDestroy().
 */
class PISMDM {
public:
#ifdef PISM_USE_TR1
  typedef std::tr1::shared_ptr<PISMDM> Ptr;
  typedef std::tr1::weak_ptr<PISMDM> WeakPtr;
#else
  typedef std::shared_ptr<PISMDM> Ptr;
  typedef std::weak_ptr<PISMDM> WeakPtr;
#endif
  PISMDM(DM dm);
  ~PISMDM();
  DM get() const;
private:
  DM m_dm;
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

  It computes grid parameters for the fine and equally-spaced vertical grid
  used in the conservation of energy and age equations.

  \section computational_grid Organization of PISM's computational grid

  PISM uses the class IceGrid to manage computational grids and their
  parameters.

  Computational grids PISM can use are
  - rectangular,
  - equally spaced in the horizintal (X and Y) directions,
  - distributed across processors in horizontal dimensions only (every column
  is stored on one processor only),
  - are periodic in both X and Y directions (in the topological sence).

  Each processor "owns" a rectangular patch of xm times ym grid points with
  indices starting at xs and ys in the X and Y directions respectively.

  The typical code performing a point-wise computation will look like

  \code
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
  for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
  // compute something at i,j
  }
  }
  \endcode

  For finite difference (and some other) computations we often need to know
  values at map-plane neighbors of a grid point. 

  We say that a patch owned by a processor is surrounded by a strip of "ghost"
  grid points belonging to patches next to the one in question. This lets us to
  access (read) values at all the eight neighbors of a grid point for \e all
  the grid points, including ones at an edge of a processor patch \e and at an
  edge of a computational domain.

  All the values \e written to ghost points will be lost next time ghost values
  are updated.

  Sometimes it is beneficial to update ghost values locally (for instance when
  a computation A uses finite differences to compute derivatives of a quantity
  produced using a purely local (point-wise) computation B). In this case the
  double loop above can be modified to look like

  \code
  for (PointsWithGhosts p(grid, ghost_width); p; p.next()) {
    const int i = p.i(), j = p.j();
    field(i,j) = value;
  }
  \endcode

  to iterate over points without ghosts, do

  \code
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    field(i,j) = value;
  }
  \endcode

*/
class IceGrid {
public:
  IceGrid(MPI_Comm c, const Config &config);
  ~IceGrid();

  PetscErrorCode report_parameters();

  PetscErrorCode allocate();  // destructor checks if DA was created, and destroys
  PetscErrorCode set_vertical_levels(const std::vector<double> &z_levels);
  PetscErrorCode compute_vertical_levels();
  PetscErrorCode compute_horizontal_spacing();
  void compute_point_neighbors(double x, double y,
                               int &i, int &j);
  std::vector<double> compute_interp_weights(double x, double y);

  void check_parameters();

  void compute_nprocs();
  void compute_ownership_ranges();
  PetscErrorCode compute_viewer_size(int target, int &x, int &y);
  PetscErrorCode printInfo(int verbosity); 
  PetscErrorCode printVertLevels(int verbosity); 
  unsigned int kBelowHeight(double height);
  PetscErrorCode create_viewer(int viewer_size, const std::string &title, PetscViewer &viewer);
  double radius(int i, int j);
  PetscErrorCode get_dm(int dm_dof, int stencil_width, PISMDM::Ptr &result);
  double convert(double, const std::string &, const std::string &) const;
  UnitSystem get_unit_system() const;

  const Config &config;
  MPI_Comm    com;
  int rank, size;
  int    xs,               //!< starting x-index of a processor sub-domain
    xm,                         //!< number of grid points (in the x-direction) in a processor sub-domain
    ys,                         //!< starting y-index of a processor sub-domain
    ym; //!< number of grid points (in the y-direction) in a processor sub-domain

  std::vector<double> zlevels; //!< vertical grid levels in the ice; correspond to the storage grid

  std::vector<double> x,             //!< x-coordinates of grid points
    y;                          //!< y-coordinates of grid points

  // Fine vertical grid and the interpolation setup:
  std::vector<double> zlevels_fine;   //!< levels of the fine vertical grid in the ice
  double   dz_fine;                    //!< spacing of the fine vertical grid
  unsigned int Mz_fine;          //!< number of levels of the fine vertical grid in the ice

  // Array ice_storage2fine contains indices of the ice storage vertical grid
  // that are just below a level of the fine grid. I.e. ice_storage2fine[k] is
  // the storage grid level just below fine-grid level k (zlevels_fine[k]).
  // Similarly for other arrays below.
  std::vector<int> ice_storage2fine, ice_fine2storage;

  SpacingType ice_vertical_spacing;
  Periodicity periodicity;
  double dzMIN,            //!< minimal vertical spacing of the storage grid in the ice
    dzMAX;                      //!< maximal vertical spacing of the storage grid in the ice

  double x0,               //!< x-coordinate of the grid center
    y0;                         //!< y-coordinate of the grid center

  double Lx, //!< half width of the ice model grid in x-direction (m)
    Ly;           //!< half width of the ice model grid in y-direction (m)

  int    Mx, //!< number of grid points in the x-direction
    My;      //!< number of grid points in the y-direction

  int    Nx, //!< number of processors in the x-direction
    Ny;      //!< number of processors in the y-direction

  std::vector<PetscInt> procs_x, //!< \brief array containing lenghts (in the x-direction) of processor sub-domains
    procs_y; //!< \brief array containing lenghts (in the y-direction) of processor sub-domains

  double dx,               //!< horizontal grid spacing
    dy;                         //!< horizontal grid spacing

  double Lz;      //!< max extent of the ice in z-direction (m)

  unsigned int Mz; //!< number of grid points in z-direction in the ice

  unsigned int initial_Mz; //!< initial number of vertical grid levels; used by the grid extension code

  Time *time;               //!< The time management object (hides calendar computations)

  //! @brief Check if a point `(i,j)` is in the strip of `stripwidth`
  //! meters around the edge of the computational domain.
  inline bool in_null_strip(int i, int j, double strip_width) {
    if (strip_width < 0.0) {
      return false;
    }
    return (x[i] <= x[0] + strip_width || x[i] >= x[Mx-1] - strip_width ||
            y[j] <= y[0] + strip_width || y[j] >= y[My-1] - strip_width);
  }
protected:
  std::map<int,PISMDM::WeakPtr> m_dms;
  double m_lambda;         //!< quadratic vertical spacing parameter
  UnitSystem m_unit_system;

  PetscErrorCode get_dzMIN_dzMAX_spacingtype();
  PetscErrorCode compute_horizontal_coordinates();
  PetscErrorCode compute_fine_vertical_grid();
  PetscErrorCode init_interpolation();

  PetscErrorCode create_dm(int da_dof, int stencil_width, DM &result);

  int dm_key(int, int);
  PetscErrorCode init_calendar(std::string &result);
private:
  // Hide copy constructor / assignment operator.
  IceGrid(IceGrid const &);
  IceGrid & operator=(IceGrid const &);
};

/** Iterator class for traversing the grid, including ghost points.
 *
 * Usage:
 *
 * `for (PointsWithGhosts p(grid, stencil_width); p; p.next()) { ... }`
 */
class PointsWithGhosts {
public:
  PointsWithGhosts(IceGrid &g, unsigned int stencil_width = 1) {
    m_i_first = g.xs - stencil_width;
    m_i_last  = g.xs + g.xm + stencil_width - 1;
    m_j_first = g.ys - stencil_width;
    m_j_last  = g.ys + g.ym + stencil_width - 1;

    m_i = m_i_first;
    m_j = m_j_first;
    m_done = false;
  }

  int i() const { return m_i; }
  int j() const { return m_j; }

  void next() {
    assert(m_done == false);
    m_j += 1;
    if (m_j > m_j_last) {
      m_j = m_j_first;        // wrap around
      m_i += 1;
    }
    if (m_i > m_i_last) {
      m_i = m_i_first;        // ensure that indexes are valid
      m_done = true;
    }
  }

  operator bool() const { return m_done == false; }
private:
  int m_i, m_j;
  int m_i_first, m_i_last, m_j_first, m_j_last;
  bool m_done;
};

/** Iterator class for traversing the grid (without ghost points).
 *
 * Usage:
 *
 * `for (Points p(grid); p; p.next()) { double foo = p.i(); ... }`
 */
class Points : public PointsWithGhosts {
public:
  Points(IceGrid &g) : PointsWithGhosts(g, 0) {}
};

} // end of namespace pism

#endif  /* __grid_hh */
