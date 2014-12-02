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

#include <petscfix.h>
#include "IceGrid.hh"
#include "pism_const.hh"
#include "PISMTime.hh"
#include "PISMTime_Calendar.hh"
#include "PISMConfig.hh"
#include "pism_options.hh"
#include <cstdlib>              // abort()
#include "error_handling.hh"

namespace pism {

PISMDM::PISMDM(DM dm) {
  m_dm = dm;
}

PISMDM::~PISMDM() {
  PetscErrorCode ierr = DMDestroy(&m_dm); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    // We can't do anything about this failure. We can't recover
    // from it, and it is almost certainly caused by a programming
    // error. So, we call abort().
    abort();
  }
}

DM PISMDM::get() const {
  return m_dm;
}

PISMDM::operator DM() const {
  return m_dm;
}

IceGrid::IceGrid(MPI_Comm c, const Config &conf)
  : config(conf), com(c), m_unit_system(config.get_unit_system()) {

  MPI_Comm_rank(com, &rank);
  MPI_Comm_size(com, &size);

  // The grid in symmetric with respect to zero by default.
  x0 = 0.0;
  y0 = 0.0;

  // initialize these data members to get rid of a valgrind warning;
  // correct values will be set in IceGrid::allocate()
  m_xs = 0;
  m_ys = 0;
  m_xm = 0;
  m_ym = 0;

  std::string word = config.get_string("grid_periodicity");
  if (word == "none") {
    periodicity = NONE;
  } else if (word == "x") {
    periodicity = X_PERIODIC;
  } else if (word == "y") {
    periodicity = Y_PERIODIC;
  } else if (word == "xy") {
    periodicity = XY_PERIODIC;
  } else {
    throw RuntimeError::formatted("grid periodicity type '%s' is invalid.",
                                  word.c_str());
  }

  word = config.get_string("grid_ice_vertical_spacing");
  if (word == "quadratic") {
    ice_vertical_spacing = QUADRATIC;
  } else if (word == "equal") {
    ice_vertical_spacing = EQUAL;
  } else {
    throw RuntimeError::formatted("ice vertical spacing type '%s' is invalid.",
                                  word.c_str());
  }

  Lx  = config.get("grid_Lx");
  Ly  = config.get("grid_Ly");
  Lz  = config.get("grid_Lz");

  m_lambda = config.get("grid_lambda");

  Mx  = static_cast<int>(config.get("grid_Mx"));
  My  = static_cast<int>(config.get("grid_My"));
  Mz  = static_cast<int>(config.get("grid_Mz"));

  Nx = Ny = 0;                  // will be set to a correct value in allocate()

  Mz_fine = 0;

  compute_vertical_levels();
  compute_horizontal_spacing();

  std::string calendar;
  PetscErrorCode ierr = init_calendar(calendar);
  if (ierr != 0) {
    throw RuntimeError("Calendar initialization failed.");
  }

  if (calendar == "360_day" || calendar == "365_day" || calendar == "noleap" || calendar == "none") {
    time = new Time(com, config, calendar, m_unit_system);
  } else {
    time = new Time_Calendar(com, config, calendar, m_unit_system);
  }
  // time->init() will be called later (in IceModel::set_grid_defaults() or
  // PIO::get_grid()).
}

/*! @brief Initialize a uniform, shallow (3 z-levels), doubly periodic grid
 * with half-widths (Lx,Ly) and Mx by My nodes.
 */
IceGrid::Ptr IceGrid::Shallow(MPI_Comm c, const Config &config,
                              double my_Lx, double my_Ly,
                              unsigned int Mx, unsigned int My, Periodicity p) {

  return IceGrid::Create(c, config, my_Lx, my_Ly, config.get("grid_Lz"),
                         Mx, My, 3, p);
}

IceGrid::Ptr IceGrid::Create(MPI_Comm c, const Config &config,
                             double my_Lx, double my_Ly, double my_Lz,
                             unsigned int Mx, unsigned int My, unsigned int Mz,
                             Periodicity p) {

  Ptr result(new IceGrid(c, config));

  result->Lx = my_Lx;
  result->Ly = my_Ly;
  result->Lz = my_Lz;
  result->Mx = Mx;
  result->My = My;
  result->periodicity = p;
  result->Mz = Mz;

  result->compute_nprocs();
  result->compute_ownership_ranges();
  result->compute_vertical_levels();
  result->compute_horizontal_spacing();
  result->allocate();

  return result;
}

IceGrid::Ptr IceGrid::Create(MPI_Comm c, const Config &config) {

  Ptr result(new IceGrid(c, config));

  // use defaults from config

  result->compute_nprocs();
  result->compute_ownership_ranges();
  result->compute_vertical_levels();
  result->compute_horizontal_spacing();
  result->allocate();

  return result;
}

/**
 * Select a calendar using the "calendar" configuration parameter, the
 * "-calendar" command-line option, or the "calendar" attribute of the
 * "time" variable in the file specified using "-time_file".
 *
 * @param[out] result selected calendar string
 *
 * @return 0 on success
 */
PetscErrorCode IceGrid::init_calendar(std::string &result) {
  PetscErrorCode ierr;

  // Set the default calendar using the config. parameter or the
  // "-calendar" option:
  result = config.get_string("calendar");

  // Check if -time_file was set and override the setting above if the
  // "calendar" attribute is found.
  std::string time_file_name;
  bool time_file_set;
  ierr = OptionsString("-time_file", "name of the file specifying the run duration",
                       time_file_name, time_file_set); CHKERRQ(ierr);
  if (time_file_set) {
    PIO nc(*this, "netcdf3");    // OK to use netcdf3

    nc.open(time_file_name, PISM_READONLY);
    {
      std::string time_name = config.get_string("time_dimension_name");
      bool time_exists = nc.inq_var(time_name);
      if (time_exists) {
        std::string tmp = nc.get_att_text(time_name, "calendar");
        if (tmp.empty() == false) {
          result = tmp;
        }
      }
    }
    nc.close();
  }

  return 0;
}

IceGrid::~IceGrid() {
  delete time;
}

//! \brief Set the vertical levels in the ice according to values in Mz, Lz,
//! and the ice_vertical_spacing data member.
/*!
Sets `dzMIN` and `dzMAX`.  Sets and re-allocates `zlevels[]`.

Uses `Mz`, `Lz`, and `ice_vertical_spacing`.  (Note that `ice_vertical_spacing`
cannot be UNKNOWN.)

This procedure is only called when a grid is determined from scratch, %e.g.
by a derived class or when bootstrapping from 2D data only, but not when
reading a model state input file (which will have its own grid,
which may not even be a grid created by this routine).
  - When `vertical_spacing` == EQUAL, the vertical grid in the ice is equally spaced:
    `zlevels[k] = k dzMIN` where `dzMIN = Lz / (Mz - 1)`.
    In this case `dzMIN = dzMAX`.
  - When `vertical_spacing` == QUADRATIC, the spacing is a quadratic function.  The intent
    is that the spacing is smaller near the base than near the top.  In particular, if
    \f$\zeta_k = k / (\mathtt{Mz} - 1)\f$ then `zlevels[k] = Lz *
    ((\f$\zeta_k\f$ / \f$\lambda\f$) * (1.0 + (\f$\lambda\f$ - 1.0)
    * \f$\zeta_k\f$))` where \f$\lambda\f$ = 4.  The value \f$\lambda\f$
    indicates the slope of the quadratic function as it leaves the base.
    Thus a value of \f$\lambda\f$ = 4 makes the spacing about four times finer
    at the base than equal spacing would be.
 */
PetscErrorCode  IceGrid::compute_vertical_levels() {

  if (Mz < 2) {
    throw RuntimeError("IceGrid::compute_ice_vertical_levels(): Mz must be at least 2.");
  }

  if (Lz <= 0) {
    throw RuntimeError("IceGrid::compute_ice_vertical_levels(): Lz must be positive.");
  }

  // Fill the levels in the ice:
  zlevels.resize(Mz);

  switch (ice_vertical_spacing) {
  case EQUAL: {
    dzMIN = Lz / ((double) Mz - 1);
    dzMAX = dzMIN;

    // Equal spacing
    for (unsigned int k=0; k < Mz - 1; k++) {
      zlevels[k] = dzMIN * ((double) k);
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    break;
  }
  case QUADRATIC: {
    // this quadratic scheme is an attempt to be less extreme in the fineness near the base.
    for (unsigned int k=0; k < Mz - 1; k++) {
      const double zeta = ((double) k) / ((double) Mz - 1);
      zlevels[k] = Lz * ((zeta / m_lambda) * (1.0 + (m_lambda - 1.0) * zeta));
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    dzMIN = zlevels[1] - zlevels[0];
    dzMAX = zlevels[Mz-1] - zlevels[Mz-2];
    break;
  }
  default:
    throw RuntimeError("IceGrid::compute_ice_vertical_levels(): ice_vertical_spacing can not be UNKNOWN.");
  }

  PetscErrorCode ierr = compute_fine_vertical_grid(); CHKERRQ(ierr);

  return 0;
}

//! Print to stdout information on computational domain and grid (other than the vertical levels themselves).
PetscErrorCode IceGrid::printInfo(const int verbosity) {
  verbPrintf(verbosity,com,
             "  IceGrid parameters:\n");
  verbPrintf(verbosity,com,
             "            Lx = %6.2f km, Ly = %6.2f km, Lz = %6.2f m,\n",
             Lx/1000.0,Ly/1000.0,Lz);
  verbPrintf(verbosity,com,
             "            x0 = %6.2f km, y0 = %6.2f km,   (coordinates of center)\n",
             x0/1000.0, y0/1000.0);
  verbPrintf(verbosity,com,
             "            Mx = %d, My = %d, Mz = %d,\n",
             Mx,My,Mz);
  verbPrintf(verbosity,com,
             "            dx = %6.3f km, dy = %6.3f km, year = %s,\n",
             dx/1000.0,dy/1000.0, time->date().c_str());
  verbPrintf(verbosity,com,
             "            Nx = %d, Ny = %d]\n",
             Nx,Ny);
  return 0;
}


//! Print the vertical levels in `zlevels[]` to stdout.
PetscErrorCode IceGrid::printVertLevels(const int verbosity) {
  verbPrintf(verbosity,com,
             "    vertical levels in ice (Mz=%d,Lz=%5.4f): ",Mz,Lz);
  for (unsigned int k=0; k < Mz; k++) {
    verbPrintf(verbosity,com," %5.4f,",zlevels[k]);
  }
  verbPrintf(verbosity,com,"\n");
  return 0;
}


//! Return the index `k` into `zlevels[]` so that `zlevels[k] <= height < zlevels[k+1]` and `k < Mz`.
unsigned int IceGrid::kBelowHeight(double height) {
  if (height < 0.0 - 1.0e-6) {
    PetscPrintf(PETSC_COMM_SELF,
                "IceGrid kBelowHeight(), rank %d, height = %5.4f is below base of ice (height must be non-negative)\n", rank, height);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  if (height > Lz + 1.0e-6) {
    PetscPrintf(PETSC_COMM_SELF,
                "IceGrid kBelowHeight(): rank %d, height = %5.4f is above top of computational grid Lz = %5.4f\n", rank, height, Lz);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  unsigned int mcurr = 0;
  while (zlevels[mcurr+1] < height) {
    mcurr++;
  }
  return mcurr;
}

//! \brief From given vertical grid zlevels[], determine `dzMIN`, `dzMAX`, and
//! determine whether ice vertical spacings are equal.
/*! The standard for equal vertical spacing in the ice is \f$10^{-8}\f$ m max
  difference between `dzMIN` and `dzMAX`.
 */
PetscErrorCode IceGrid::get_dzMIN_dzMAX_spacingtype() {

  // ice:
  dzMIN = Lz;
  dzMAX = 0.0;
  for (unsigned int k = 0; k < Mz - 1; k++) {
    const double mydz = zlevels[k+1] - zlevels[k];
    dzMIN = PetscMin(mydz,dzMIN);
    dzMAX = PetscMax(mydz,dzMAX);
  }
  if (PetscAbs(dzMAX - dzMIN) <= 1.0e-8) {
    ice_vertical_spacing = EQUAL;
  } else {
    ice_vertical_spacing = UNKNOWN;
  }

  return 0;
}

//! \brief Computes the number of processors in the X- and Y-directions.
void IceGrid::compute_nprocs() {

  if (My <= 0) {
    throw RuntimeError("'My' is invalid.");
  }

  Nx = (int)(0.5 + sqrt(((double)Mx)*((double)size)/((double)My)));

  if (Nx == 0) {
    Nx = 1;
  }

  while (Nx > 0) {
    Ny = size/Nx;
    if (Nx*Ny == size) {
      break;
    }
    Nx--;
  }

  if (Mx > My && Nx < Ny) {int _Nx = Nx; Nx = Ny; Ny = _Nx;}

  if ((Mx / Nx) < 2) {          // note: integer division
    throw RuntimeError::formatted("Can't distribute a %d x %d grid across %d processors!",
                                  Mx, My, size);
  }

  if ((My / Ny) < 2) {          // note: integer division
    throw RuntimeError::formatted("Can't distribute a %d x %d grid across %d processors!",
                                  Mx, My, size);
  }
}

//! \brief Computes processor ownership ranges corresponding to equal area
//! distribution among processors.
/*!
 * Expects grid.Nx and grid.Ny to be valid.
 */
void IceGrid::compute_ownership_ranges() {

  procs_x.resize(Nx);
  procs_y.resize(Ny);

  for (int i=0; i < Nx; i++) {
    procs_x[i] = Mx/Nx + ((Mx % Nx) > i);
  }

  for (int i=0; i < Ny; i++) {
    procs_y[i] = My/Ny + ((My % Ny) > i);
  }
}

//! \brief Create the PETSc DM for the horizontal grid. Determine how
//! the horizontal grid is divided among processors.
/*!
  This procedure should only be called after the parameters describing the
  horizontal computational box (Lx,Ly) and the parameters for the horizontal
  grid (Mx,My) are already determined. In particular, the input file (either \c
  -i or `-boot_file`) and user options (like `-Mx`) must have already been
  read to determine the parameters, and any conflicts must have been resolved.

  This method contains the "fundamental" transpose: "My,Mx" instead of "Mx,My"
  in the DMDACreate2d call; this transpose allows us to index arrays by "[i][j]"
  (where 'i' corresponds to 'x' and 'j' to 'y') and be consistent about
  meanings of 'x', 'y', 'u' and 'v'.

  Unfortunately this means that PETSc viewers appear transposed.

  This choice should be virtually invisible, unless you're using DALocalInfo
  structures.

  \note PETSc order: x in columns, y in rows, indexing as array[y][x]. PISM
  order: x in rows, y in columns, indexing as array[x][y].
 */
PetscErrorCode IceGrid::allocate() {

  check_parameters();

  unsigned int max_stencil_width = (int)config.get("grid_max_stencil_width");

  try {
    PISMDM::Ptr tmp = this->get_dm(1, max_stencil_width);
  } catch (RuntimeError) {
    throw RuntimeError::formatted("can't distribute the %d x %d grid across %d processors.",
                                  Mx, My, size);
  }

  // hold on to a DM corresponding to dof=1, stencil_width=0 (it will
  // be needed for I/O operations)
  m_dm_scalar_global = this->get_dm(1, 0);

  DMDALocalInfo info;
  DMDAGetLocalInfo(*m_dm_scalar_global, &info);
  // this continues the fundamental transpose
  m_xs = info.ys; m_xm = info.ym;
  m_ys = info.xs; m_ym = info.xm;

  return 0;
}

//! Sets grid vertical levels; sets Mz and Lz from input.  Checks input for consistency.
PetscErrorCode IceGrid::set_vertical_levels(const std::vector<double> &new_zlevels) {

  if (new_zlevels.size() < 2) {
    throw RuntimeError("IceGrid::set_vertical_levels(): Mz has to be at least 2.");
  }

  if ((!is_increasing(new_zlevels)) || (PetscAbs(new_zlevels[0]) > 1.0e-10)) {
    throw RuntimeError("IceGrid::set_vertical_levels(): invalid zlevels; must be strictly increasing and start with z=0.");
  }

  Mz  =  (int)new_zlevels.size();
  Lz  =  new_zlevels.back();

  zlevels  = new_zlevels;

  get_dzMIN_dzMAX_spacingtype();
  compute_fine_vertical_grid();

  return 0;
}


//! Compute horizontal spacing parameters `dx` and `dy` using `Mx`, `My`, `Lx`, `Ly` and periodicity.
/*!
The grid used in PISM, in particular the PETSc DAs used here, are periodic in x and y.
This means that the ghosted values ` foo[i+1][j], foo[i-1][j], foo[i][j+1], foo[i][j-1]`
for all 2D Vecs, and similarly in the x and y directions for 3D Vecs, are always available.
That is, they are available even if i,j is a point at the edge of the grid.  On the other
hand, by default, `dx`  is the full width  `2 * Lx`  divided by  `Mx - 1`.
This means that we conceive of the computational domain as starting at the `i = 0`
grid location and ending at the  `i = Mx - 1`  grid location, in particular.
This idea is not quite compatible with the periodic nature of the grid.

The upshot is that if one computes in a truly periodic way then the gap between the
`i = 0`  and  `i = Mx - 1`  grid points should \em also have width  `dx`.
Thus we compute  `dx = 2 * Lx / Mx`.
 */
PetscErrorCode IceGrid::compute_horizontal_spacing() {

  if (periodicity & X_PERIODIC) {
    dx = 2.0 * Lx / Mx;
  } else {
    dx = 2.0 * Lx / (Mx - 1);
  }

  if (periodicity & Y_PERIODIC) {
    dy = 2.0 * Ly / My;
  } else {
    dy = 2.0 * Ly / (My - 1);
  }

  compute_horizontal_coordinates();

  return 0;
}

//! Computes fine vertical spacing in the ice.
/*! The computations in IceModel::temperatureStep(), IceModel::enthalpyDrainageStep()
  and IceModel::ageStep() use a fine equally-spaced grid.

  This method computes the number of levels and the levels themselves so that
  fine equally-spaced grids in ice and bedrock have the same spacing (if
  bedrock is present).

  Mapping to and from the storage grid occurs in IceModelVec3 methods.

  Note that the computational grid in the ice is allowed to exceed Lz; we need
  this to match spacings (see above). The temperature field is extrapolated to
  the extra level using the value from the topmost level.
 */
PetscErrorCode IceGrid::compute_fine_vertical_grid() {

  // the smallest of the spacings used in ice and bedrock:
  double my_dz_fine = dzMIN;

  Mz_fine = static_cast<int>(ceil(Lz / my_dz_fine) + 1);
  my_dz_fine = Lz / (Mz_fine - 1);

  // both ice and bedrock will have this spacing
  dz_fine = my_dz_fine;

  // allocate arrays:
  zlevels_fine.resize(Mz_fine);

  // compute levels in the ice:
  for (unsigned int k = 0; k < Mz_fine; k++) {
    zlevels_fine[k] = ((double) k) * dz_fine;
  }
  // Note that it's allowed to go over Lz.

  init_interpolation();

  return 0;
}

//! Fills arrays ice_storage2fine, ice_fine2storage, bed_storage2fine,
//! bed_fine2storage with indices of levels that are just below
PetscErrorCode IceGrid::init_interpolation() {
  unsigned int m = 0;

  // ice: storage -> fine
  ice_storage2fine.resize(Mz_fine);
  m = 0;
  for (unsigned int k = 0; k < Mz_fine; k++) {
    if (zlevels_fine[k] >= Lz) {
      ice_storage2fine[k] = Mz - 1;
      continue;
    }

    while (zlevels[m + 1] < zlevels_fine[k]) {
      m++;
    }

    ice_storage2fine[k] = m;
  }

  // ice: fine -> storage
  ice_fine2storage.resize(Mz);
  m = 0;
  for (unsigned int k = 0; k < Mz; k++) {
    while (m < Mz_fine - 1 && zlevels_fine[m + 1] < zlevels[k]) {
      m++;
    }

    ice_fine2storage[k] = m;
  }

  return 0;
}

//! \brief Computes values of x and y corresponding to the computational grid,
//! with accounting for periodicity.
PetscErrorCode IceGrid::compute_horizontal_coordinates() {

  x.resize(Mx);
  y.resize(My);

  // Here x_min, x_max define the extent of the computational domain,
  // which is not necessarily the same thing as the smallest and
  // largest values of x.
  double
    x_min = x0 - Lx,
    x_max = x0 + Lx;
  if (periodicity & X_PERIODIC) {
    for (int i = 0; i < Mx; ++i) {
      x[i] = x_min + (i + 0.5) * dx;
    }
    x[Mx - 1] = x_max - 0.5*dx;
  } else {
    for (int i = 0; i < Mx; ++i) {
      x[i] = x_min + i * dx;
    }
    x[Mx - 1] = x_max;
  }

  double
    y_min = y0 - Ly,
    y_max = y0 + Ly;
  if (periodicity & Y_PERIODIC) {
    for (int i = 0; i < My; ++i) {
      y[i] = y_min + (i + 0.5) * dy;
    }
    y[My - 1] = y_max - 0.5*dy;
  } else {
    for (int i = 0; i < My; ++i) {
      y[i] = y_min + i * dy;
    }
    y[My - 1] = y_max;
  }

  return 0;
}

//! \brief Returns the distance from the point (i,j) to the origin.
double IceGrid::radius(int i, int j) {
  return sqrt(PetscSqr(x[i]) + PetscSqr(y[j]));
}

//! \brief Report grid parameters.
PetscErrorCode IceGrid::report_parameters() {

  verbPrintf(2,com, "computational domain and grid:\n");

  // report on grid
  verbPrintf(2,com,
             "                grid size   %d x %d x %d\n",
             Mx,My,Mz);
  // report on computational box
  verbPrintf(2,com,
             "           spatial domain   %.2f km x %.2f km x %.2f m\n",
             2*Lx/1000.0,2*Ly/1000.0,Lz);

  // report on grid cell dims
  verbPrintf(2,com,
             "     horizontal grid cell   %.2f km x %.2f km\n",
             dx/1000.0,dy/1000.0);
  if (ice_vertical_spacing == EQUAL) {
    verbPrintf(2,com,
               "  vertical spacing in ice   dz = %.3f m (equal spacing)\n",
               dzMIN);
  } else {
    verbPrintf(2,com,
               "  vertical spacing in ice   uneven, %d levels, %.3f m < dz < %.3f m\n",
               Mz, dzMIN, dzMAX);
  }
  verbPrintf(3,com,
             "   fine spacing used in energy/age   fMz = %d, fdz = %.3f m\n",
             Mz_fine, dz_fine);
  if (Mz_fine > 1000) {
    verbPrintf(2,com,
               "\n\nWARNING: Using more than 1000 ice vertical levels internally in energy/age computation!\n\n");
  }

  // report on time axis
  //   FIXME:  this could use pism_config:summary_time_unit_name instead of fixed "years"
  verbPrintf(2, com,
             "   time interval (length)   [%s, %s]  (%s years, using the '%s' calendar)\n",
             time->start_date().c_str(),
             time->end_date().c_str(),
             time->run_length().c_str(),
             time->calendar().c_str());

  // if -verbose (=-verbose 3) then (somewhat redundantly) list parameters of grid
  printInfo(3);

  verbPrintf(5,com,
             "  REALLY verbose output on IceGrid:\n");
  printVertLevels(5);

  return 0;
}

//! Computes the size of a diagnostic viewer window.
PetscErrorCode IceGrid::compute_viewer_size(int target_size, int &X, int &Y) {

  // aim for smaller dimension equal to target, larger dimension larger by Ly/Lx or Lx/Ly proportion
  const double  yTOx = Ly / Lx;
  if (Ly > Lx) {
    X = target_size;
    Y = (int) ((double)target_size * yTOx);
  } else {
    Y = target_size;
    X = (int) ((double)target_size / yTOx);
  }

  // if either dimension is larger than twice the target, shrink appropriately
  if (X > 2 * target_size) {
    Y = (int) ((double)(Y) * (2.0 * (double)target_size / (double)(X)));
    X = 2 * target_size;
  } else if (Y > 2 * target_size) {
    X = (int) ((double)(X) * (2.0 * (double)target_size / (double)(Y)));
    Y = 2 * target_size;
  }

  // make sure minimum dimension is sufficient to see
  if (X < 20) {
    X = 20;
  }
  if (Y < 20) {
    Y = 20;
  }
  return 0;
}

//! Creates a run-time diagnostic viewer.
PetscErrorCode IceGrid::create_viewer(int viewer_size, const std::string &title, PetscViewer &viewer) {
  PetscErrorCode ierr;
  int X, Y;

  compute_viewer_size(viewer_size, X, Y);

  // note we reverse x <-> y; see IceGrid::allocate() for original reversal
  ierr = PetscViewerDrawOpen(com, NULL, title.c_str(),
                             PETSC_DECIDE, PETSC_DECIDE, Y, X, &viewer);
  PISM_PETSC_CHK(ierr, "PetscViewerDrawOpen");

  // following should be redundant, but may put up a title even under 2.3.3-p1:3 where
  // there is a no titles bug
  PetscDraw draw;
  ierr = PetscViewerDrawGetDraw(viewer, 0, &draw);
  PISM_PETSC_CHK(ierr, "PetscViewerDrawGetDraw");
  ierr = PetscDrawSetTitle(draw, title.c_str());
  PISM_PETSC_CHK(ierr, "PetscDrawSetTitle");

  return 0;
}

//! \brief Computes indices of grid points to the lower left and upper right from (X,Y).
/*!
 * \code
 * 3       2
 * o-------o
 * |       |
 * |    +  |
 * o-------o
 * 0       1
 * \endcode
 *
 * If "+" is the point (X,Y), then (i_left, j_bottom) corresponds to
 * point "0" and (i_right, j_top) corresponds to point "2".
 *
 * Does not check if the resulting indexes are in the current
 * processor's domain. Ensures that computed indexes are within the
 * grid.
 */
void IceGrid::compute_point_neighbors(double X, double Y,
                                      int &i_left, int &i_right,
                                      int &j_bottom, int &j_top) {
  i_left = (int)floor((X - x[0])/dx);
  j_bottom = (int)floor((Y - y[0])/dy);

  i_right = i_left + 1;
  j_top = j_bottom + 1;

  if (i_left < 0) {
    i_left = i_right;
  }

  if (i_right > Mx - 1) {
    i_right = i_left;
  }

  if (j_bottom < 0) {
    j_bottom = j_top;
  }

  if (j_top > My - 1) {
    j_top = j_bottom;
  }
}

//! \brief Compute 4 interpolation weights necessary for linear interpolation
//! from the current grid. See compute_point_neighbors for the ordering of
//! neighbors.
std::vector<double> IceGrid::compute_interp_weights(double X, double Y) {
  int i_left = 0, i_right = 0, j_bottom = 0, j_top = 0;
  // these values (zeros) are used when interpolation is impossible
  double alpha = 0.0, beta = 0.0;

  compute_point_neighbors(X, Y, i_left, i_right, j_bottom, j_top);

  if (i_left != i_right) {
    assert(x[i_right] - x[i_left] != 0.0);
    alpha = (X - x[i_left]) / (x[i_right] - x[i_left]);
  }

  if (j_bottom != j_top) {
    assert(y[j_top] - y[j_bottom] != 0.0);
    beta  = (Y - x[j_bottom]) / (y[j_top] - y[j_bottom]);
  }

  std::vector<double> result(4);
  result[0] = alpha * beta;
  result[1] = (1 - alpha) * beta;
  result[2] = (1 - alpha) * (1 - beta);
  result[3] = alpha * (1 - beta);

  return result;
}

//! \brief Checks grid parameters usually set at bootstrapping for validity.
void IceGrid::check_parameters() {

  if (Mx < 3) {
    throw RuntimeError("Mx has to be at least 3.");
  }

  if (My < 3) {
    throw RuntimeError("My has to be at least 3.");
  }

  if (Mz < 2) {
    throw RuntimeError("Mz must be at least 2.");
  }

  if (Lx <= 0) {
    throw RuntimeError("Lx has to be positive.");
  }

  if (Ly <= 0) {
    throw RuntimeError("Ly has to be positive.");
  }

  if (Lz <= 0) {
    throw RuntimeError("Lz must be positive.");
  }

  // A single record of a time-dependent variable cannot exceed 2^32-4
  // bytes in size. See the NetCDF User's Guide
  // <http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#g_t64-bit-Offset-Limitations>.
  // Here we use "long int" to avoid integer overflow.
  const long int two_to_thirty_two = 4294967296L;
  const long int Mx_long = Mx, My_long = My, Mz_long = Mz;
  if (Mx_long * My_long * Mz_long * sizeof(double) > two_to_thirty_two - 4 &&
      ((config.get_string("output_format") == "netcdf3") ||
       (config.get_string("output_format") == "pnetcdf"))) {
    throw RuntimeError::formatted("The computational grid is too big to fit in a NetCDF-3 file.\n"
                                  "Each 3D variable requires %lu Mb.\n"
                                  "Please use '-o_format quilt' or re-build PISM with parallel NetCDF-4 or HDF5\n"
                                  "and use '-o_format netcdf4_parallel' or '-o_format hdf5' to proceed.",
                                  Mx_long * My_long * Mz_long * sizeof(double) / (1024 * 1024));
  }

}

PISMDM::Ptr IceGrid::get_dm(int da_dof, int stencil_width) {
  PISMDM::Ptr result;

  if (da_dof < 0 || da_dof > 10000) {
    throw RuntimeError::formatted("Invalid da_dof argument: %d", da_dof);
  }

  if (stencil_width < 0 || stencil_width > 10000) {
    throw RuntimeError::formatted("Invalid stencil_width argument: %d", stencil_width);
  }

  int j = this->dm_key(da_dof, stencil_width);

  if (m_dms[j].expired() == true) {
    DM tmp = this->create_dm(da_dof, stencil_width);

    result = PISMDM::Ptr(new PISMDM(tmp));
    m_dms[j] = result;
  } else {
    result = m_dms[j].lock();
  }

  return result;
}

UnitSystem IceGrid::get_unit_system() const {
  return m_unit_system;
}

double IceGrid::convert(double value, const std::string &unit1, const std::string &unit2) const {
  return m_unit_system.convert(value, unit1, unit2);
}

DM IceGrid::create_dm(int da_dof, int stencil_width) {
  DM result;

  verbPrintf(3, com,
             "* Creating a DM with dof=%d and stencil_width=%d...\n",
             da_dof, stencil_width);

  PetscErrorCode ierr;
  ierr = DMDACreate2d(com,
#if PETSC_VERSION_LT(3,5,0)
                      DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
#else
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
#endif
                      DMDA_STENCIL_BOX,
                      My, Mx, // N, M
                      Ny, Nx, // n, m
                      da_dof, stencil_width,
                      &procs_y[0], &procs_x[0], // ly, lx
                      &result);
  PISM_PETSC_CHK(ierr,"DMDACreate2d");

  return result;
}

// Computes the key corresponding to the DM with given dof and stencil_width.
int IceGrid::dm_key(int da_dof, int stencil_width) {
  return 10000 * da_dof + stencil_width;
}

int IceGrid::xs() const {
  return m_xs;
}

int IceGrid::ys() const {
  return m_ys;
}

int IceGrid::xm() const {
  return m_xm;
}

int IceGrid::ym() const {
  return m_ym;
}

} // end of namespace pism
