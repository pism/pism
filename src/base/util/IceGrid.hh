// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __grid_hh
#define __grid_hh

#include <petscdmda.h>
#include <vector>
#include <string>

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

class PISMTime;
class PISMProf;
class NCConfigVariable;

typedef enum {UNKNOWN = 0, EQUAL, QUADRATIC} SpacingType;
typedef enum {NONE = 0, NOT_PERIODIC =0, X_PERIODIC = 1, Y_PERIODIC = 2, XY_PERIODIC = 3} Periodicity;

//! Describes the PISM grid and the distribution of data across processors.
/*!
  This class holds parameters describing the grid, including the vertical
  spacing and which part of the horizontal grid is owned by the processor. It
  contains the dimensions of the PISM (4-dimensional, x*y*z*time) computational
  box. The vertical spacing can be quite arbitrary.

  It creates and destroys a two dimensional \c PETSc \c DA (distributed array).
  The creation of this \c DA is the point at which PISM gets distributed across
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
  int GHOSTS = 1;
  for (int i=grid.xs - GHOSTS; i<grid.xs+grid.xm + GHOSTS; ++i) {
    for (int j=grid.ys - GHOSTS; j<grid.ys+grid.ym + GHOSTS; ++j) {
    // compute something at i,j
    }
  }
  \endcode

 */
class IceGrid {
public:
  IceGrid(MPI_Comm c, PetscMPIInt r, PetscMPIInt s, const NCConfigVariable &config);
  ~IceGrid();

  PetscErrorCode report_parameters();

  PetscErrorCode createDA();  // destructor checks if DA was created, and destroys
  PetscErrorCode createDA(int procs_x, int procs_y,
			  int* &lx, int* &ly);
  PetscErrorCode set_vertical_levels(vector<double> z_levels);
  PetscErrorCode compute_vertical_levels();
  PetscErrorCode compute_horizontal_spacing();
  void compute_point_neighbors(PetscReal x, PetscReal y,
                               int &i, int &j);
  vector<PetscReal> compute_interp_weights(PetscReal x, PetscReal y);

  void check_parameters();

  void compute_nprocs();
  void compute_ownership_ranges();
  PetscErrorCode compute_viewer_size(int target, int &x, int &y);
  PetscErrorCode printInfo(int verbosity); 
  PetscErrorCode printVertLevels(int verbosity); 
  int       kBelowHeight(PetscScalar height);
  PetscErrorCode create_viewer(int viewer_size, string title, PetscViewer &viewer);
  PetscReal      radius(int i, int j);

  const NCConfigVariable &config;
  MPI_Comm    com;
  PetscMPIInt rank, size;
  DM          da2;

  int    xs,               //!< starting x-index of a processor sub-domain
    xm,                         //!< number of grid points (in the x-direction) in a processor sub-domain
    ys,                         //!< starting y-index of a processor sub-domain
    ym; //!< number of grid points (in the y-direction) in a processor sub-domain

  vector<double> zlevels; //!< vertical grid levels in the ice; correspond to the storage grid

  vector<double> x,             //!< x-coordinates of grid points
    y;          //!< y-coordinates of grid points

  // Fine vertical grid and the interpolation setup:
  vector<double> zlevels_fine;   //!< levels of the fine vertical grid in the ice
  PetscReal   dz_fine;                    //!< spacing of the fine vertical grid
  int    Mz_fine;          //!< number of levels of the fine vertical grid in the ice

  // Array ice_storage2fine contains indices of the ice storage vertical grid
  // that are just below a level of the fine grid. I.e. ice_storage2fine[k] is
  // the storage grid level just below fine-grid level k (zlevels_fine[k]).
  // Similarly for other arrays below.
  vector<int> ice_storage2fine, ice_fine2storage;

  SpacingType ice_vertical_spacing;
  Periodicity periodicity;
  PetscScalar dzMIN,            //!< minimal vertical spacing of the storage grid in the ice
    dzMAX;                      //!< maximal vertical spacing of the storage grid in the ice

  PetscScalar x0,               //!< x-coordinate of the grid center
    y0;	   //!< y-coordinate of the grid center

  PetscScalar Lx, //!< half width of the ice model grid in x-direction (m)
    Ly; //!< half width of the ice model grid in y-direction (m)

  int    Mx, //!< number of grid points in the x-direction
    My; //!< number of grid points in the y-direction

  int    Nx, //!< number of processors in the x-direction
    Ny; //!< number of processors in the y-direction

  vector<int> procs_x, //!< \brief array containing lenghts (in the x-direction) of processor sub-domains
    procs_y; //!< \brief array containing lenghts (in the y-direction) of processor sub-domains

  PetscScalar dx,               //!< horizontal grid spacing
    dy;                         //!< horizontal grid spacing

  PetscScalar Lz;      //!< extent of the ice in z-direction (m)

  int    Mz; //!< number of grid points in z-direction in the ice

  int initial_Mz; //!< initial number of vertical grid levels; used by the grid extension code

  int max_stencil_width;   //!< \brief maximum stencil width supported by
                                //!< the DA in this IceGrid object

  PISMProf *profiler;           //!< PISM profiler object; allows tracking how long a computation takes
  PISMTime *time;               //!< The time management object (hides calendar computations)
protected:
  PetscScalar lambda;	 //!< quadratic vertical spacing parameter
  PetscErrorCode get_dzMIN_dzMAX_spacingtype();
  PetscErrorCode compute_horizontal_coordinates();
  PetscErrorCode compute_fine_vertical_grid();
  PetscErrorCode init_interpolation();

private:
  // Hide copy constructor / assignment operator.
  IceGrid(IceGrid const &);
  IceGrid & operator=(IceGrid const &);

};

#endif	/* __grid_hh */

