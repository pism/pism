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

#include <petscfix.h>
#include "IceGrid.hh"
#include "pism_const.hh"
#include "PISMTime.hh"
#include "PISMGregorianTime.hh"
#include "PISMProf.hh"
#include "NCVariable.hh"


IceGrid::IceGrid(MPI_Comm c, PetscMPIInt r, PetscMPIInt s,
		 const NCConfigVariable &conf)
  : config(conf), com(c), rank(r), size(s){ 

  // The grid in symmetric with respect to zero by default.
  x0 = 0.0;
  y0 = 0.0;

  string word = config.get_string("grid_periodicity");
  if (word == "none")
    periodicity = NONE;
  else if (word == "x")
    periodicity = X_PERIODIC;
  else if (word == "y")
    periodicity = Y_PERIODIC;
  else if (word == "xy")
    periodicity = XY_PERIODIC;
  else {
    PetscPrintf(com, 
		"ERROR: grid periodicity type '%s' is invalid.\n",
		word.c_str());
    PISMEnd();
  }

  word = config.get_string("grid_ice_vertical_spacing");
  if (word == "quadratic")
    ice_vertical_spacing = QUADRATIC;
  else if (word == "equal")
    ice_vertical_spacing = EQUAL;
  else {
    PetscPrintf(com, 
		"ERROR: ice vertical spacing type '%s' is invalid.\n",
		word.c_str());
    PISMEnd();
  }

  Lx  = config.get("grid_Lx");
  Ly  = config.get("grid_Ly");
  Lz  = config.get("grid_Lz");

  lambda = config.get("grid_lambda");

  Mx  = static_cast<PetscInt>(config.get("grid_Mx"));
  My  = static_cast<PetscInt>(config.get("grid_My"));
  Mz  = static_cast<PetscInt>(config.get("grid_Mz"));

  Nx = Ny = 0;			// will be set to a correct value in createDA()
  initial_Mz = 0;		// will be set to a correct value in
				// IceModel::check_maximum_thickness()

  max_stencil_width = 2;
  da2 = PETSC_NULL;

  Mz_fine = 0;

  compute_vertical_levels();
  compute_horizontal_spacing();

  profiler = new PISMProf(com, rank, size);

  if (config.get_string("calendar") == "gregorian") {
    time = new PISMGregorianTime(com, config);
  } else {
    time = new PISMTime(com, config);
  }
  // time->init() will be called later (in IceModel::set_grid_defaults() or
  // PIO::get_grid()).
}


IceGrid::~IceGrid() {
  if (da2 != PETSC_NULL) {
    DMDestroy(&da2);
  }

  delete time;
  delete profiler;
}

//! \brief Set the vertical levels in the ice according to values in Mz, Lz,
//! and the ice_vertical_spacing data member.
/*!
Sets \c dzMIN and \c dzMAX.  Sets and re-allocates \c zlevels[].

Uses \c Mz, \c Lz, and \c ice_vertical_spacing.  (Note that \c ice_vertical_spacing
cannot be UNKNOWN.)

This procedure is only called when a grid is determined from scratch, %e.g. 
by a derived class or when bootstrapping from 2D data only, but not when 
reading a model state input file (which will have its own grid,
which may not even be a grid created by this routine).
  - When \c vertical_spacing == EQUAL, the vertical grid in the ice is equally spaced:
    <tt>zlevels[k] = k dzMIN</tt> where <tt>dzMIN = Lz / (Mz - 1)</tt>.  
    In this case <tt>dzMIN = dzMAX</tt>.
  - When \c vertical_spacing == QUADRATIC, the spacing is a quadratic function.  The intent 
    is that the spacing is smaller near the base than near the top.  In particular, if 
    \f$\zeta_k = k / (\mathtt{Mz} - 1)\f$ then <tt>zlevels[k] = Lz * 
    ( (\f$\zeta_k\f$ / \f$\lambda\f$) * (1.0 + (\f$\lambda\f$ - 1.0) 
    * \f$\zeta_k\f$) )</tt> where \f$\lambda\f$ = 4.  The value \f$\lambda\f$ 
    indicates the slope of the quadratic function as it leaves the base.  
    Thus a value of \f$\lambda\f$ = 4 makes the spacing about four times finer 
    at the base than equal spacing would be.
 */
PetscErrorCode  IceGrid::compute_vertical_levels() {
  
  if (Mz < 2) {
    SETERRQ(com, 2,"IceGrid::compute_ice_vertical_levels(): Mz must be at least 2.");
  }

  if (Lz <= 0) {
    SETERRQ(com, 4, "IceGrid::compute_ice_vertical_levels(): Lz must be positive.");
  }

  // Fill the levels in the ice:
  zlevels.resize(Mz);

  switch (ice_vertical_spacing) {
  case EQUAL: {
    dzMIN = Lz / ((PetscScalar) Mz - 1);
    dzMAX = dzMIN;
    
    // Equal spacing
    for (PetscInt k=0; k < Mz - 1; k++) {
      zlevels[k] = dzMIN * ((PetscScalar) k);
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    break;
  }
  case QUADRATIC: {
    // this quadratic scheme is an attempt to be less extreme in the fineness near the base.
    for (PetscInt k=0; k < Mz - 1; k++) {
      const PetscScalar zeta = ((PetscScalar) k) / ((PetscScalar) Mz - 1);
      zlevels[k] = Lz * ( (zeta / lambda) * (1.0 + (lambda - 1.0) * zeta) );
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    dzMIN = zlevels[1] - zlevels[0];
    dzMAX = zlevels[Mz-1] - zlevels[Mz-2];
    break;
  }
  default:
    SETERRQ(com, 1,"IceGrid::compute_ice_vertical_levels(): ice_vertical_spacing can not be UNKNOWN.");
  }

  PetscErrorCode ierr = compute_fine_vertical_grid(); CHKERRQ(ierr);

  return 0;
}

//! Print to stdout information on computational domain and grid (other than the vertical levels themselves).
PetscErrorCode IceGrid::printInfo(const int verbosity) {
  PetscErrorCode ierr;
  ierr = verbPrintf(verbosity,com,
         "  IceGrid parameters:\n"); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            Lx = %6.2f km, Ly = %6.2f km, Lz = %6.2f m,\n",
         Lx/1000.0,Ly/1000.0,Lz); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            x0 = %6.2f km, y0 = %6.2f km,   (coordinates of center)\n",
		    x0/1000.0, y0/1000.0); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            Mx = %d, My = %d, Mz = %d,\n",
         Mx,My,Mz); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            dx = %6.3f km, dy = %6.3f km, year = %s,\n",
                    dx/1000.0,dy/1000.0, time->date().c_str()); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            Nx = %d, Ny = %d]\n",
         Nx,Ny); CHKERRQ(ierr);
  return 0;
}


//! Print the vertical levels in \c zlevels[] to stdout.
PetscErrorCode IceGrid::printVertLevels(const int verbosity) {
  PetscErrorCode ierr;
  ierr = verbPrintf(verbosity,com,
     "    vertical levels in ice (Mz=%d,Lz=%5.4f): ",Mz,Lz); CHKERRQ(ierr);
  for (PetscInt k=0; k < Mz; k++) {
    ierr = verbPrintf(verbosity,com," %5.4f,",zlevels[k]); CHKERRQ(ierr);
  }
  ierr = verbPrintf(verbosity,com,"\n"); CHKERRQ(ierr);
  return 0;
}


//! Return the index \c k into \c zlevels[] so that <tt>zlevels[k] <= height < zlevels[k+1]</tt> and <tt>k < Mz</tt>.
PetscInt IceGrid::kBelowHeight(PetscScalar height) {
  if (height < 0.0 - 1.0e-6) {
    PetscPrintf(com, 
       "IceGrid kBelowHeight(): height = %5.4f is below base of ice (height must be non-negative)\n",
       height);
    PISMEnd();
  }
  if (height > Lz + 1.0e-6) {
    PetscPrintf(com, 
       "IceGrid kBelowHeight(): height = %5.4f is above top of computational grid Lz = %5.4f\n",
       height,Lz);
    PISMEnd();
  }
  PetscInt mcurr = 0;
//  while ((zlevels[mcurr+1] <= height) && (mcurr+1 < Mz)) {
  while (zlevels[mcurr+1] < height) {
    mcurr++;
  }
  return mcurr;
}

//! \brief From given vertical grid zlevels[], determine \c dzMIN, \c dzMAX, and
//! determine whether ice vertical spacings are equal.
/*! The standard for equal vertical spacing in the ice is \f$10^{-8}\f$ m max
  difference between \c dzMIN and \c dzMAX.
 */
PetscErrorCode IceGrid::get_dzMIN_dzMAX_spacingtype() {

  // ice:
  dzMIN = Lz; 
  dzMAX = 0.0;
  for (PetscInt k = 0; k < Mz - 1; k++) {
    const PetscScalar mydz = zlevels[k+1] - zlevels[k];
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
    PetscPrintf(com, "PISM ERROR: 'My' is invalid. Exiting...\n");
    PISMEnd();
  }

  Nx = (PetscInt)(0.5 + sqrt(((double)Mx)*((double)size)/((double)My)));

  if (Nx == 0) Nx = 1;

  while (Nx > 0) {
    Ny = size/Nx;
    if (Nx*Ny == size) break;
    Nx--;
  }

  if (Mx > My && Nx < Ny) {PetscInt _Nx = Nx; Nx = Ny; Ny = _Nx;}

  if ((Mx / Nx) < 2) {		// note: integer division
    PetscPrintf(com, "PISM ERROR: Can't distribute a %d x %d grid across %d processors!\n",
		Mx, My, size);
    PISMEnd();
  }

  if ((My / Ny) < 2) {		// note: integer division
    PetscPrintf(com, "PISM ERROR: Can't distribute a %d x %d grid across %d processors!\n",
		Mx, My, size);
    PISMEnd();
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
  
  for (PetscInt i=0; i < Nx; i++) {
    procs_x[i] = Mx/Nx + ((Mx % Nx) > i);
  }

  for (PetscInt i=0; i < Ny; i++) {
    procs_y[i] = My/Ny + ((My % Ny) > i);
  }
}

//! \brief Create the PETSc DA \c da2 for the horizontal grid. Determine how
//! the horizontal grid is divided among processors.
/*!
  This procedure should only be called after the parameters describing the
  horizontal computational box (Lx,Ly) and the parameters for the horizontal
  grid (Mx,My) are already determined. In particular, the input file (either \c
  -i or \c -boot_file) and user options (like \c -Mx) must have already been
  read to determine the parameters, and any conflicts must have been resolved.

  This method contains the "fundamental" transpose: "My,Mx" instead of "Mx,My"
  in the DACreate2d call; this transpose allows us to index arrays by "[i][j]"
  (where 'i' corresponds to 'x' and 'j' to 'y') and be consistent about
  meanings of 'x', 'y', 'u' and 'v'.

  Unfortunately this means that PETSc viewers appear transposed.

  This choice should be virtually invisible, unless you're using DALocalInfo
  structures.

  See IceModelVec3 for creation of three-dimensional DAs.

  \note PETSc order: x in columns, y in rows, indexing as array[y][x]. PISM
  order: x in rows, y in columns, indexing as array[x][y].
 */
PetscErrorCode IceGrid::createDA() {
  PetscErrorCode ierr;

  if (Mx < 3) {
    SETERRQ(com, 1, "IceGrid::createDA(): Mx has to be at least 3.");
  }

  if (My < 3) {
    SETERRQ(com, 2, "IceGrid::createDA(): My has to be at least 3.");
  }

  if (Lx <= 0) {
    SETERRQ(com, 3, "IceGrid::createDA(): Lx has to be positive.");
  }

  if (Ly <= 0) {
    SETERRQ(com, 3, "IceGrid::createDA(): Ly has to be positive.");
  }

  if (da2 != PETSC_NULL)
    SETERRQ(com, 1, "IceGrid::createDA(): da2 != PETSC_NULL");

  // Transpose:
  ierr = DMDACreate2d(com,
                      DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      My, Mx,
                      Ny, Nx,
                      1, max_stencil_width, // dof, stencil width
                      &procs_y[0], &procs_x[0],
                      &da2);
  if (ierr != 0) {
    PetscErrorCode ierr2;
    ierr2 = verbPrintf(1, com,
		       "PISM ERROR: can't distribute the %d x %d grid across %d processors...\n"
		       "Exiting...\n", Mx, My, ierr);
    CHKERRQ(ierr2);
    PISMEnd();
  }

  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da2, &info); CHKERRQ(ierr);
  // this continues the fundamental transpose
  xs = info.ys; xm = info.ym;
  ys = info.xs; ym = info.xm;

  profiler->Nx = Nx;
  profiler->Ny = Ny;
  profiler->set_grid_size(xm*ym);

  return 0;
}

//! \brief A version allowing to specify the grid distribution across processors.
/*!
 * \param my_procs_x Number of processors in the X-direction
 * \param my_procs_y Number of processors in the Y-direction
 * \param lx The array containing numbers of grid points per processor in the X-direction
 * \param ly The array containing numbers of grid points per processor in the Y-direction
 *
 * \note my_procs_x times my_procs_y has to be equal to size
 * \note lx has to add up to Mx.
 * \note ly has to add up to My.
 */
PetscErrorCode IceGrid::createDA(PetscInt my_procs_x, PetscInt my_procs_y,
				 PetscInt* &lx, PetscInt* &ly) {
  PetscErrorCode ierr;

  if (da2 != PETSC_NULL) {
    ierr = DMDestroy(&da2); CHKERRQ(ierr);
  }

  // Transpose:
  ierr = DMDACreate2d(com,
                      DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      My, Mx,
                      my_procs_y, my_procs_x,
                      1, max_stencil_width, // dof, stencil width
                      ly, lx, &da2);
  if (ierr != 0) {
    PetscErrorCode ierr2;
    ierr2 = verbPrintf(1, com,
		       "PISM ERROR: can't distribute the %d x %d grid across %d processors...\n"
		       "Exiting...\n", Mx, My, size);
    CHKERRQ(ierr2);
    PISMEnd();
  }
  
 CHKERRQ(ierr);

  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da2, &info); CHKERRQ(ierr);
  // this continues the fundamental transpose
  xs = info.ys; xm = info.ym;
  ys = info.xs; ym = info.xm;

  profiler->Nx = Nx;
  profiler->Ny = Ny;

  return 0;
}

//! Sets grid vertical levels; sets Mz and Lz from input.  Checks input for consistency.
PetscErrorCode IceGrid::set_vertical_levels(vector<double> new_zlevels) {
  PetscErrorCode ierr;

  if (new_zlevels.size() < 2) {
    SETERRQ(com, 1, "IceGrid::set_vertical_levels(): Mz has to be at least 2.");
  }

  if ( (!is_increasing(new_zlevels)) || (PetscAbs(new_zlevels[0]) > 1.0e-10) ) {
    SETERRQ(com, 3, "IceGrid::set_vertical_levels(): invalid zlevels; must be strictly increasing and start with z=0.");
  }

  Mz  =  (int)new_zlevels.size();
  Lz  =  new_zlevels.back();

  zlevels  = new_zlevels;

  get_dzMIN_dzMAX_spacingtype();
  ierr = compute_fine_vertical_grid(); CHKERRQ(ierr);

  return 0;
}


//! Compute horizontal spacing parameters \c dx and \c dy using \c Mx, \c My, \c Lx, \c Ly and periodicity.
/*! 
The grid used in PISM, in particular the PETSc DAs used here, are periodic in x and y.
This means that the ghosted values <tt> foo[i+1][j], foo[i-1][j], foo[i][j+1], foo[i][j-1]</tt>
for all 2D Vecs, and similarly in the x and y directions for 3D Vecs, are always available.
That is, they are available even if i,j is a point at the edge of the grid.  On the other
hand, by default, \c dx  is the full width  <tt>2 * Lx</tt>  divided by  <tt>Mx - 1</tt>.
This means that we conceive of the computational domain as starting at the <tt>i = 0</tt>
grid location and ending at the  <tt>i = Mx - 1</tt>  grid location, in particular.  
This idea is not quite compatible with the periodic nature of the grid.

The upshot is that if one computes in a truly periodic way then the gap between the  
<tt>i = 0</tt>  and  <tt>i = Mx - 1</tt>  grid points should \em also have width  \c dx.  
Thus we compute  <tt>dx = 2 * Lx / Mx</tt>.
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
  PetscErrorCode ierr;

  // the smallest of the spacings used in ice and bedrock:
  PetscScalar my_dz_fine = dzMIN;

  Mz_fine = static_cast<PetscInt>(ceil(Lz / my_dz_fine) + 1);
  my_dz_fine = Lz / (Mz_fine - 1);

  // both ice and bedrock will have this spacing
  dz_fine = my_dz_fine;

  // allocate arrays:
  zlevels_fine.resize(Mz_fine);

  // compute levels in the ice:
  for (PetscInt k = 0; k < Mz_fine; k++)
    zlevels_fine[k] = ((PetscScalar) k) * dz_fine;
  // Note that it's allowed to go over Lz.

  ierr = init_interpolation(); CHKERRQ(ierr);

  return 0;
}

//! Fills arrays ice_storage2fine, ice_fine2storage, bed_storage2fine,
//! bed_fine2storage with indices of levels that are just below 
PetscErrorCode IceGrid::init_interpolation() {
  PetscInt m;

  // ice: storage -> fine
  ice_storage2fine.resize(Mz_fine);
  m = 0;
  for (PetscInt k = 0; k < Mz_fine; k++) {
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
  for (PetscInt k = 0; k < Mz; k++) {
    while (zlevels_fine[m + 1] < zlevels[k]) {
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

  double x_min = -Lx + x0,
    x_max = Lx + x0,
    y_min = -Ly + y0,
    y_max = Ly + y0;

  if (periodicity & X_PERIODIC) {
    x_max -= dx;
  }

  if (periodicity & Y_PERIODIC) {
    y_max -= dy;
  }

  for (int i = 0; i < Mx; ++i)
    x[i] = x_min + i * dx;
  x[Mx - 1] = x_max;

  for (int i = 0; i < My; ++i)
    y[i] = y_min + i * dy;
  y[My - 1] = y_max;

  return 0;
}

//! \brief Returns the distance from the point (i,j) to the origin.
PetscReal IceGrid::radius(PetscInt i, PetscInt j) {
  return sqrt(PetscSqr(x[i]) + PetscSqr(y[j]));
}

//! \brief Report grid parameters.
PetscErrorCode IceGrid::report_parameters() {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,com, "computational domain and grid:\n"); CHKERRQ(ierr);

  // report on grid
  ierr = verbPrintf(2,com, 
           "                grid size   %d x %d x %d\n",
           Mx,My,Mz); CHKERRQ(ierr);
  // report on computational box
  ierr = verbPrintf(2,com, 
           "           spatial domain   %.2f km x %.2f km x %.2f m\n",
           2*Lx/1000.0,2*Ly/1000.0,Lz); CHKERRQ(ierr);
  
  // report on grid cell dims
  ierr = verbPrintf(2,com, 
           "     horizontal grid cell   %.2f km x %.2f km\n",
                    dx/1000.0,dy/1000.0); CHKERRQ(ierr);
  if (ice_vertical_spacing == EQUAL) {
    ierr = verbPrintf(2,com, 
           "  vertical spacing in ice   dz = %.3f m (equal spacing)\n",
                    dzMIN); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,com, 
           "  vertical spacing in ice   uneven, %d levels, %.3f m < dz < %.3f m\n",
		    Mz, dzMIN, dzMAX); CHKERRQ(ierr);
  }
  ierr = verbPrintf(3,com, 
                    "   fine spacing used in energy/age   fMz = %d, fdz = %.3f m\n",
                    Mz_fine, dz_fine); CHKERRQ(ierr);
  if (Mz_fine > 1000) {
    ierr = verbPrintf(2,com,
      "\n\nWARNING: Using more than 1000 ice vertical levels internally in energy/age computation!\n\n");
      CHKERRQ(ierr);
  }

  // report on time axis
  ierr = verbPrintf(2, com,
           "   time interval (length)   [%s, %s]  (%s years)\n",
		    time->start_date().c_str(),
                    time->end_date().c_str(),
                    time->run_length().c_str()); CHKERRQ(ierr);

  // if -verbose (=-verbose 3) then (somewhat redundantly) list parameters of grid
  ierr = printInfo(3); CHKERRQ(ierr);

  ierr = verbPrintf(5,com,
       "  REALLY verbose output on IceGrid:\n"); CHKERRQ(ierr);
  ierr = printVertLevels(5); CHKERRQ(ierr);
  
  return 0;
}

//! Computes the size of a diagnostic viewer window.
PetscErrorCode IceGrid::compute_viewer_size(int target_size, int &X, int &Y) {

  // aim for smaller dimension equal to target, larger dimension larger by Ly/Lx or Lx/Ly proportion
  const double  yTOx = Ly / Lx;
  if (Ly > Lx) {
    X = target_size; 
    Y = (PetscInt) ((double)target_size * yTOx); 
  } else {
    Y = target_size; 
    X = (PetscInt) ((double)target_size / yTOx);
  }
  
  // if either dimension is larger than twice the target, shrink appropriately
  if (X > 2 * target_size) {
    Y = (PetscInt) ( (double)(Y) * (2.0 * (double)target_size / (double)(X)) );
    X = 2 * target_size;
  } else if (Y > 2 * target_size) {
    X = (PetscInt) ( (double)(X) * (2.0 * (double)target_size / (double)(Y)) );
    Y = 2 * target_size;
  }
  
  // make sure minimum dimension is sufficient to see
  if (X < 20)   X = 20;
  if (Y < 20)   Y = 20;
  return 0;
}

//! Creates a run-time diagnostic viewer.
PetscErrorCode IceGrid::create_viewer(PetscInt viewer_size, string title, PetscViewer &viewer) {
  PetscErrorCode ierr;
  int X, Y;

  ierr = compute_viewer_size(viewer_size, X, Y); CHKERRQ(ierr);

  // note we reverse x <-> y; see IceGrid::createDA() for original reversal
  ierr = PetscViewerDrawOpen(com, PETSC_NULL, title.c_str(),
			     PETSC_DECIDE, PETSC_DECIDE, Y, X, &viewer);  CHKERRQ(ierr);

  // following should be redundant, but may put up a title even under 2.3.3-p1:3 where
  // there is a no titles bug
  PetscDraw draw;
  ierr = PetscViewerDrawGetDraw(viewer, 0, &draw); CHKERRQ(ierr);
  ierr = PetscDrawSetTitle(draw, title.c_str()); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes indices of a grid point to the lower left of the point (x,y).
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
 * If "+" is the point (X,Y), this method returns indices of the grid point
 * "0".
 *
 * Does not if the resulting i and j are valid or in the current processor's
 * domain.
 */
void IceGrid::compute_point_neighbors(PetscReal X, PetscReal Y,
                                      PetscInt &i, PetscInt &j) {
  i = (int)floor((X - x[0])/dx);
  j = (int)floor((Y - y[0])/dy);
}

//! \brief Compute 4 interpolation weights necessary for linear interpolation
//! from the current grid. See compute_point_neighbors for the ordering of
//! neighbors.
vector<PetscReal> IceGrid::compute_interp_weights(PetscReal X, PetscReal Y) {
  vector<PetscReal> result;
  int i,j;
  double alpha, beta;

  compute_point_neighbors(X, Y, i, j);

  result.resize(4);

  if ((i >= 0) && (i < Mx))
    alpha = (X - x[i]) / dx;
  else alpha = 0;

  if ((i >= 0) && (j < My))
    beta  = (Y - x[j]) / dy;
  else beta = 0;

  result[0] = alpha * beta;
  result[1] = (1 - alpha) * beta;
  result[2] = (1 - alpha) * (1 - beta);
  result[3] = alpha * (1 - beta);

  return result;
}

//! \brief Checks grid parameters usually set at bootstrapping for validity.
void IceGrid::check_parameters() {

  if (Mx < 3) {
    PetscPrintf(com, "PISM ERROR: Mx has to be at least 3.\n");
    PISMEnd();
  }

  if (My < 3) {
    PetscPrintf(com, "PISM ERROR: My has to be at least 3.\n");
    PISMEnd();
  }

  if (Mz < 2) {
    PetscPrintf(com, "PISM ERROR: Mz must be at least 2.\n");
    PISMEnd();
  }

  if (Lx <= 0) {
    PetscPrintf(com, "PISM ERROR: Lx has to be positive.\n");
    PISMEnd();
  }

  if (Ly <= 0) {
    PetscPrintf(com, "PISM ERROR: Ly has to be positive.\n");
    PISMEnd();
  }

  if (Lz <= 0) {
    PetscPrintf(com, "PISM ERROR: Lz must be positive.\n");
    PISMEnd();
  }

  if (Mx * My * Mz * sizeof(double) > pow(2, 32) - 4 &&
      config.get_string("output_format") != "netcdf4_parallel") {
    PetscPrintf(com,
                "PISM ERROR: The computational grid is too big to fit in a NetCDF-3 file.\n"
                "            Each 3D variable requires %lu Mb.\n"
                "            Please re-build PISM with parallel NetCDF-4\n"
                "            and use \"-o_format netcdf4_parallel\" to proceed.\n",
                Mx * My * Mz * sizeof(double) / (1024 * 1024));
    PISMEnd();
  }

}

