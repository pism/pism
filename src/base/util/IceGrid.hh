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

#include <cassert>

#include "pism_const.hh"

namespace pism {

class Time;
class Prof;
class Config;

typedef enum {UNKNOWN = 0, EQUAL, QUADRATIC} SpacingType;
typedef enum {NONE = 0, NOT_PERIODIC = 0, X_PERIODIC = 1, Y_PERIODIC = 2, XY_PERIODIC = 3} Periodicity;

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
  operator DM() const;
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

#ifdef PISM_USE_TR1
  typedef std::tr1::shared_ptr<IceGrid> Ptr;
#else
  typedef std::shared_ptr<IceGrid> Ptr;
#endif

  static Ptr Shallow(MPI_Comm c, const Config &config,
                     double my_Lx, double my_Ly,
                     unsigned int Mx, unsigned int My, Periodicity p);

  static Ptr Create(MPI_Comm c, const Config &config,
                    double my_Lx, double my_Ly, double my_Lz,
                    unsigned int Mx, unsigned int My, unsigned int Mz,
                    Periodicity p);

  static Ptr Create(MPI_Comm c, const Config &config);

  // static Ptr Bootstrapping(MPI_Comm c, const Config &config,
  //                          const std::string &filename);

  void report_parameters() const;

  // only of these two should be called:
  void set_vertical_levels(const std::vector<double> &z_levels);
  void compute_vertical_levels();

  void allocate();  // FIXME! allocate in the constructor!

  void compute_point_neighbors(double X, double Y,
                               int &i_left, int &i_right,
                               int &j_bottom, int &j_top);
  std::vector<double> compute_interp_weights(double x, double y);

  unsigned int kBelowHeight(double height);
  PISMDM::Ptr get_dm(int dm_dof, int stencil_width);

  // FIXME: these should be moved into a "Context" class
  double convert(double, const std::string &, const std::string &) const;
  UnitSystem get_unit_system() const;
  double radius(int i, int j);

  //! Starting x-index of a processor sub-domain
  int xs() const;
  //! Number of grid points (in the x-direction) in a processor sub-domain
  int xm() const;
  //! Starting y-index of a processor sub-domain
  int ys() const;
  //! Number of grid points (in the y-direction) in a processor sub-domain
  int ym() const;

  const std::vector<double>& x() const;
  double x(size_t i) const;

  const std::vector<double>& y() const;
  double y(size_t i) const;

  const std::vector<double>& z() const;
  double z(size_t i) const;

  double dx() const;
  double dy() const;

  unsigned int Mx() const;
  unsigned int My() const;
  unsigned int Mz() const;

  double Lx() const;
  double Ly() const;
  double Lz() const;
  // FIXME: remove this
  void set_Lx(double Lx) {m_Lx = Lx;}
  void set_Ly(double Ly) {m_Ly = Ly;}
  void set_Lz(double Lz) {m_Lz = Lz;}

  double x0() const;
  double y0() const;
  // FIXME: remove this
  void set_x0(double x0) {m_x0 = x0;}
  void set_y0(double y0) {m_y0 = y0;}

  // FIXME: remove this
  void set_Mx(unsigned int Mx) {m_Mx = Mx;}
  void set_My(unsigned int My) {m_My = My;}
  void set_Mz(unsigned int Mz) {m_Mz = Mz;}

  double dz_min() const;
  double dz_max() const;

  Profiling profiling;

  const Config &config;
  const MPI_Comm com;
  // int to match types used by MPI
  int rank, size;

  // Fine vertical grid and the interpolation setup:
  std::vector<double> zlevels_fine;   //!< levels of the fine vertical grid in the ice
  double dz_fine;                    //!< spacing of the fine vertical grid
  unsigned int Mz_fine;          //!< number of levels of the fine vertical grid in the ice

  // Array ice_storage2fine contains indices of the ice storage vertical grid
  // that are just below a level of the fine grid. I.e. ice_storage2fine[k] is
  // the storage grid level just below fine-grid level k (zlevels_fine[k]).
  // Similarly for other arrays below.
  std::vector<int> ice_storage2fine, ice_fine2storage;

  SpacingType ice_vertical_spacing;
  Periodicity periodicity;

  // FIXME: this should be moved into a "Context" class
  Time *time;               //!< The time management object (hides calendar computations)
private:
  unsigned int m_Nx, //!< number of processors in the x-direction
    m_Ny;      //!< number of processors in the y-direction

  //! @brief array containing lenghts (in the x-direction) of processor sub-domains
  std::vector<PetscInt> m_procs_x;
  //! @brief array containing lenghts (in the y-direction) of processor sub-domains
  std::vector<PetscInt> m_procs_y;

  std::vector<double> m_x,             //!< x-coordinates of grid points
    m_y;                          //!< y-coordinates of grid points

  //! vertical grid levels in the ice; correspond to the storage grid
  std::vector<double> m_zlevels;

  int m_xs, m_xm, m_ys, m_ym;
  double m_dx,               //!< horizontal grid spacing
    m_dy;                    //!< horizontal grid spacing
  //! number of grid points in the x-direction
  unsigned int m_Mx;
  //! number of grid points in the y-direction
  unsigned int m_My;
  //! number of grid points in z-direction in the ice
  unsigned int m_Mz;

  //! x-coordinate of the grid center
  double m_x0;
  //! y-coordinate of the grid center
  double m_y0;

  //! half width of the ice model grid in x-direction (m)
  double m_Lx;
  //! half width of the ice model grid in y-direction (m)
  double m_Ly;
  //! max extent of the ice in z-direction (m)
  double m_Lz;

  //! minimal vertical spacing of the storage grid in the ice
  double m_dz_min;
  //! maximal vertical spacing of the storage grid in the ice
  double m_dz_max;

  std::map<int,PISMDM::WeakPtr> m_dms;
  double m_lambda;         //!< quadratic vertical spacing parameter
  const UnitSystem m_unit_system;

  // This DM is used for I/O operations and is not owned by any
  // IceModelVec (so far, anyway). We keep a pointer to it here to
  // avoid re-allocating it many times.
  PISMDM::Ptr m_dm_scalar_global;

  void check_parameters();

  void ownership_ranges_from_options();
  void compute_nprocs();
  void compute_ownership_ranges();

  void get_dz_min_dz_max_spacingtype();
  void compute_horizontal_spacing();
  void compute_horizontal_coordinates();
  void compute_fine_vertical_grid();
  void init_interpolation();

  DM create_dm(int da_dof, int stencil_width);

  int dm_key(int, int);
  std::string init_calendar();

  // Hide copy constructor / assignment operator.
  IceGrid(IceGrid const &);
  IceGrid & operator=(IceGrid const &);
};

//! @brief Check if a point `(i,j)` is in the strip of `stripwidth`
//! meters around the edge of the computational domain.
inline bool in_null_strip(const IceGrid& grid, int i, int j, double strip_width) {
  return (strip_width >= 0.0                               &&
          (grid.x(i)  <= grid.x(0) + strip_width           ||
           grid.x(i)  >= grid.x(grid.Mx()-1) - strip_width ||
           grid.y(j)  <= grid.y(0) + strip_width           ||
           grid.y(j)  >= grid.y(grid.My()-1) - strip_width));
}

/** Iterator class for traversing the grid, including ghost points.
 *
 * Usage:
 *
 * `for (PointsWithGhosts p(grid, stencil_width); p; p.next()) { ... }`
 */
class PointsWithGhosts {
public:
  PointsWithGhosts(IceGrid &g, unsigned int stencil_width = 1) {
    m_i_first = g.xs() - stencil_width;
    m_i_last  = g.xs() + g.xm() + stencil_width - 1;
    m_j_first = g.ys() - stencil_width;
    m_j_last  = g.ys() + g.ym() + stencil_width - 1;

    m_i = m_i_first;
    m_j = m_j_first;
    m_done = false;
  }

  int i() const {
    return m_i;
  }
  int j() const {
    return m_j;
  }

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

  operator bool() const {
    return m_done == false;
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
 * `for (Points p(grid); p; p.next()) { double foo = p.i(); ... }`
 */
class Points : public PointsWithGhosts {
public:
  Points(IceGrid &g) : PointsWithGhosts(g, 0) {}
};

} // end of namespace pism

#endif  /* __grid_hh */
