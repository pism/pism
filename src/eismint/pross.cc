// Copyright (C) 2007--2011 Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include "SSAFD.hh"
#include "PISMIO.hh"
#include "Timeseries.hh"

static char help[] =
  "Driver for EISMINT-Ross diagnostic velocity computation in ice shelf.\n"
  "Illustrates use of SSA stress balance plus I/O plus time series,\n"
  "without the time-stepping mass continuity and conservation of energy\n"
  "components of PISM.\n\n";

PetscErrorCode read_riggs_and_compare(IceGrid &grid, PISMVars &vars, IceModelVec2V &vel_ssa) {
  PetscErrorCode  ierr;
  PetscTruth      riggsSet;
  char            riggsfile[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-riggs", riggsfile,
                               PETSC_MAX_PATH_LEN, &riggsSet); CHKERRQ(ierr);
  if (riggsSet == PETSC_FALSE)
    return 0;

  IceModelVec2S *longitude, *latitude;
  IceModelVec2Mask *mask;

  longitude = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
  if (longitude == NULL) SETERRQ(1, "longitude is not available");

  latitude = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (latitude == NULL) SETERRQ(1, "latitude is not available");

  mask = dynamic_cast<IceModelVec2Mask*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(1, "mask is not available");

  ierr = verbPrintf(2,grid.com,"comparing to RIGGS data in %s ...\n",
                    riggsfile); CHKERRQ(ierr);

  Timeseries latdata(&grid, "riggslat", "count"),
    londata(&grid, "riggslon", "count"),
    magdata(&grid, "riggsmag", "count"),
    udata(&grid, "riggsu", "count"),
    vdata(&grid, "riggsv", "count");
  PetscInt    len;
  PetscScalar **clat, **clon;

  ierr = longitude->get_array(clon); CHKERRQ(ierr);
  ierr =  latitude->get_array(clat); CHKERRQ(ierr);
  ierr =    vel_ssa.begin_access(); CHKERRQ(ierr);
  ierr =      mask->begin_access();  CHKERRQ(ierr);

  ierr = magdata.set_units("m year-1", ""); CHKERRQ(ierr);
  ierr =   udata.set_units("m year-1", ""); CHKERRQ(ierr);
  ierr =   vdata.set_units("m year-1", ""); CHKERRQ(ierr);

  ierr = latdata.read(riggsfile); CHKERRQ(ierr);
  ierr = londata.read(riggsfile); CHKERRQ(ierr);
  ierr = magdata.read(riggsfile); CHKERRQ(ierr);
  ierr =   udata.read(riggsfile); CHKERRQ(ierr);
  ierr =   vdata.read(riggsfile); CHKERRQ(ierr);
      
  // same length for all vars here
  len = latdata.length();
  PetscScalar  goodptcount = 0.0, ChiSqr = 0.0;
  for (PetscInt k = 0; k<len; k++) {
    PetscScalar lat, lon, mag, u, v;
    lat = latdata[k];
    lon = londata[k];
    mag = magdata[k];
    u   = udata[k];
    v   = vdata[k];
    ierr = verbPrintf(4,grid.com,
                      " RIGGS[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
                      k,lat,lon,mag,u,v); CHKERRQ(ierr); 
    const PetscScalar origdlat = (-5.42445 - (-12.3325)) / 110.0;
    const PetscScalar lowlat = -12.3325 - origdlat * 46.0;
    const PetscScalar dlat = (-5.42445 - lowlat) / (float) (grid.My - 1);        
    const PetscScalar lowlon = -5.26168;
    const PetscScalar dlon = (3.72207 - lowlon) / (float) (grid.Mx - 1);
    const int         cj = (int) floor((lat - lowlat) / dlat);
    const int         ci = (int) floor((lon - lowlon) / dlon);
    if ((ci >= grid.xs) && (ci < grid.xs+grid.xm) && (cj >= grid.ys) && (cj < grid.ys+grid.ym)) {
      const PetscScalar cu = secpera * vel_ssa(ci,cj).u;
      const PetscScalar cv = secpera * vel_ssa(ci,cj).v;
      const PetscScalar cmag = sqrt(PetscSqr(cu)+PetscSqr(cv));
      ierr = verbPrintf(4,PETSC_COMM_SELF,
                        " PISM%d[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
                        grid.rank,k,clat[ci][cj],clon[ci][cj],cmag,cu,cv); CHKERRQ(ierr); 
      if (mask->value(ci,cj) == MASK_FLOATING) {
        goodptcount += 1.0;
        ChiSqr += PetscSqr(u-cu)+PetscSqr(v-cv);
      }
    }
  }
  ChiSqr = ChiSqr / PetscSqr(30.0); // see page 48 of MacAyeal et al
  PetscScalar g_goodptcount, g_ChiSqr;
  ierr = PetscGlobalSum(&goodptcount, &g_goodptcount, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&ChiSqr, &g_ChiSqr, grid.com); CHKERRQ(ierr);
  ierr = verbPrintf(4,grid.com,"number of RIGGS data points = %d\n"
                    "number of RIGGS points in computed ice shelf region = %8.2f\n",
                    len, g_goodptcount); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"Chi^2 statistic for computed results compared to RIGGS is %10.3f\n",
                    g_ChiSqr * (156.0 / g_goodptcount)); CHKERRQ(ierr);

  ierr = longitude->end_access(); CHKERRQ(ierr);
  ierr =  latitude->end_access(); CHKERRQ(ierr);
  ierr =      mask->end_access(); CHKERRQ(ierr);
  ierr =    vel_ssa.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode compute_errors(IceGrid &grid, PISMVars &vars, IceModelVec2V &vel_ssa) {
  // average over grid, where observed velocity is accurate *according to
  // 111by147grid.dat*, the difference between computed and observed u,v
  PetscErrorCode  ierr;
  PetscScalar  uerr=0.0, verr=0.0, relvecerr=0.0, accN=0.0, 
    accArea=0.0, maxcComputed=0.0, vecErrAcc = 0.0;
  PetscScalar  **azi, **mag, **acc, **H;

  IceModelVec2S *thickness, *obsAzimuth, *obsMagnitude, *obsAccurate;
  IceModelVec2Mask *mask;

  mask = dynamic_cast<IceModelVec2Mask*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(1, "mask is not available");

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  obsAzimuth = dynamic_cast<IceModelVec2S*>(vars.get("azi_obs"));
  if (obsAzimuth == NULL) SETERRQ(1, "azi_obs is not available");

  obsMagnitude = dynamic_cast<IceModelVec2S*>(vars.get("mag_obs"));
  if (obsMagnitude == NULL) SETERRQ(1, "mag_obs is not available");

  obsAccurate = dynamic_cast<IceModelVec2S*>(vars.get("accur"));
  if (obsAccurate == NULL) SETERRQ(1, "accur is not available");

  const PetscScalar area = grid.dx * grid.dy;
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = thickness->get_array(H); CHKERRQ(ierr);
  ierr = obsAzimuth->get_array(azi); CHKERRQ(ierr);    
  ierr = obsMagnitude->get_array(mag); CHKERRQ(ierr);    
  ierr = obsAccurate->get_array(acc); CHKERRQ(ierr);    
  ierr = vel_ssa.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask->is_floating(i,j) && (H[i][j] > 1.0)) {
        const PetscScalar ccomputed = vel_ssa(i,j).magnitude();
        maxcComputed = PetscMax(maxcComputed,ccomputed);
        if (PetscAbs(acc[i][j] - 1.0) < 0.1) {
          accN += 1.0;
          accArea += area;
          const PetscScalar uobs = mag[i][j] * sin((pi/180.0) * azi[i][j]);
          const PetscScalar vobs = mag[i][j] * cos((pi/180.0) * azi[i][j]);
          // compare from readme.txt
          // uxbar(i0,j0) = magvel(i0,j0)* sin(3.1415926/180.*azvel(i0,j0))
          // uybar(i0,j0) = magvel(i0,j0)* cos(3.1415926/180.*azvel(i0,j0))
          // debug:
          //verbPrintf(1,grid.com,"i,j=%d,%d:  observed = (%5.4f,%5.4f),   computed =  (%5.4f,%5.4f)\n",
          //           i,j,uobs*secpera,vobs*secpera,vel_bar(i,j).u*secpera,vel_bar(i,j).v*secpera);
          const PetscScalar Dv = PetscAbs(vobs - vel_ssa(i,j).v);
          const PetscScalar Du = PetscAbs(uobs - vel_ssa(i,j).u);
          verr += Dv;
          uerr += Du;
          relvecerr += (PetscSqr(Dv) + PetscSqr(Du)) / (PetscSqr(vobs) + PetscSqr(uobs));
          vecErrAcc += (PetscSqr(Dv) + PetscSqr(Du)) * area;
        }
      }
    }
  }
  ierr = obsMagnitude->end_access(); CHKERRQ(ierr);
  ierr = obsAzimuth->end_access(); CHKERRQ(ierr);
  ierr = obsAccurate->end_access(); CHKERRQ(ierr);
  ierr = vel_ssa.end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  PetscScalar  guerr, gverr, grelvecerr, gaccN, 
    gaccArea, gmaxcComputed, gvecErrAcc;
  ierr = PetscGlobalMax(&maxcComputed, &gmaxcComputed, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&accN, &gaccN, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&accArea, &gaccArea, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&verr, &gverr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&uerr, &guerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&relvecerr, &grelvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&vecErrAcc, &gvecErrAcc, grid.com); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"maximum computed speed in ice shelf is %10.3f (m/a)\n",
                    gmaxcComputed * secpera); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"ERRORS relative to observations of Ross Ice Shelf:\n");
  CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  [number of grid points in 'accurate observed area' = %d]\n",
                    (int) gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  [area of 'accurate observed area' = %9.4f (km^2)]\n",
                    gaccArea / 1e6); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  following are average errors computed over 'accurate observed area':\n");
  CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  average error in x-comp of vel       = %9.3f (m/a)\n",
                    (gverr * secpera) / gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  average error in y-comp of vel       = %9.3f (m/a)\n",
                    (guerr * secpera) / gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  average relative error in vector vel = %9.5f\n",
                    grelvecerr / gaccN); CHKERRQ(ierr);
  gvecErrAcc = secpera * sqrt(gvecErrAcc) / sqrt(gaccArea);
  ierr = verbPrintf(2,grid.com,
                    "  rms average error in vector vel      = %9.3f (m/a)\n",
                    gvecErrAcc); CHKERRQ(ierr);
  return 0;
}

//! \brief Allocate necessary variables.
PetscErrorCode allocate_vars(IceGrid &grid, PISMVars &vars) {
  PetscErrorCode ierr;
  const PetscInt WIDE_STENCIL = 2;

  IceModelVec2S *obsAzimuth, *obsMagnitude, *obsAccurate,
    *thickness, *surface, *bed, *tauc, *longitude, *latitude;
  IceModelVec2Mask *mask, *bc_mask;
  IceModelVec2V *vel_bc;
  IceModelVec3 *enthalpy;

  obsAzimuth = new IceModelVec2S;
  obsMagnitude = new IceModelVec2S;
  obsAccurate = new IceModelVec2S;
  thickness = new IceModelVec2S;
  surface = new IceModelVec2S;
  bed = new IceModelVec2S;
  tauc = new IceModelVec2S;
  longitude = new IceModelVec2S;
  latitude = new IceModelVec2S;

  mask = new IceModelVec2Mask;
  bc_mask = new IceModelVec2Mask;

  vel_bc = new IceModelVec2V;

  enthalpy = new IceModelVec3;

  ierr = obsAzimuth->create(grid, "azi_obs", true); CHKERRQ(ierr);
  ierr = obsAzimuth->set_attrs("", "observed ice velocity azimuth",
                               "degrees_east", ""); CHKERRQ(ierr);
  ierr = vars.add(*obsAzimuth); CHKERRQ(ierr);

  ierr = obsMagnitude->create(grid, "mag_obs", true); CHKERRQ(ierr);
  ierr = obsMagnitude->set_attrs("", "observed ice velocity magnitude",
                                 "m s-1", ""); CHKERRQ(ierr);
  ierr = obsMagnitude->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  obsMagnitude->write_in_glaciological_units = true;
  ierr = vars.add(*obsMagnitude); CHKERRQ(ierr);

  ierr = obsAccurate->create(grid, "accur", true); CHKERRQ(ierr);
  ierr = obsAccurate->set_attrs("", "flag for accurate observed velocity",
                                "", ""); CHKERRQ(ierr);
  ierr = vars.add(*obsAccurate); CHKERRQ(ierr);

  // ice upper surface elevation
  ierr = surface->create(grid, "usurf", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = surface->set_attrs("diagnostic", "ice upper surface elevation",
                            "m", "surface_altitude"); CHKERRQ(ierr);
  ierr = vars.add(*surface); CHKERRQ(ierr);

  // land ice thickness
  ierr = thickness->create(grid, "thk", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = thickness->set_attrs("model_state", "land ice thickness",
                              "m", "land_ice_thickness"); CHKERRQ(ierr);
  ierr = thickness->set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = vars.add(*thickness); CHKERRQ(ierr);

  // bedrock surface elevation
  ierr = bed->create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = bed->set_attrs("model_state", "bedrock surface elevation",
                        "m", "bedrock_altitude"); CHKERRQ(ierr);
  ierr = vars.add(*bed); CHKERRQ(ierr);

  // yield stress for basal till (plastic or pseudo-plastic model)
  ierr = tauc->create(grid, "tauc", false); CHKERRQ(ierr);
  ierr = tauc->set_attrs("diagnostic", 
                         "yield stress for basal till (plastic or pseudo-plastic model)",
                         "Pa", ""); CHKERRQ(ierr);
  ierr = vars.add(*tauc); CHKERRQ(ierr);

  ierr = enthalpy->create(grid, "enthalpy", true); CHKERRQ(ierr);
  ierr = enthalpy->set_attrs("model_state",
                             "ice enthalpy (includes sensible heat, latent heat, pressure)",
                             "J kg-1", ""); CHKERRQ(ierr);
  ierr = vars.add(*enthalpy); CHKERRQ(ierr);

  // grounded_dragging_floating integer mask
  ierr = mask->create(grid, "mask", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = mask->set_attrs("model_state", "grounded_dragging_floating integer mask",
                         "", ""); CHKERRQ(ierr);
  vector<double> mask_values(6);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_SHEET;
  mask_values[2] = MASK_DRAGGING_SHEET;
  mask_values[3] = MASK_FLOATING;
  mask_values[4] = MASK_ICE_FREE_OCEAN;
  mask_values[5] = MASK_OCEAN_AT_TIME_0;
  ierr = mask->set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = mask->set_attr("flag_meanings",
                        "ice_free_bedrock sheet dragging_sheet floating ice_free_ocean ocean_at_time_zero");
  CHKERRQ(ierr);
  mask->output_data_type = NC_BYTE;
  ierr = vars.add(*mask); CHKERRQ(ierr);

  ierr = bc_mask->create(grid, "bcflag", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = bc_mask->set_attrs("model_state", "boundary condition locations",
                            "", ""); CHKERRQ(ierr);
  vector<double> bc_mask_values(2);
  bc_mask_values[0] = MASK_ICE_FREE_BEDROCK;
  bc_mask_values[1] = MASK_SHEET;
  ierr = bc_mask->set_attr("flag_values", bc_mask_values); CHKERRQ(ierr);
  ierr = bc_mask->set_attr("flag_meanings", "no_data bc_condition"); CHKERRQ(ierr);
  bc_mask->output_data_type = NC_BYTE;
  ierr = vars.add(*bc_mask); CHKERRQ(ierr);

  ierr = vel_bc->create(grid, "bar", false); CHKERRQ(ierr); // ubar and vbar
  ierr = vel_bc->set_attrs("intent",
                           "X-component of the SSA velocity boundary conditions",
                           "m s-1", "", 0); CHKERRQ(ierr);
  ierr = vel_bc->set_attrs("intent",
                           "Y-component of the SSA velocity boundary conditions",
                           "m s-1", "", 1); CHKERRQ(ierr);
  ierr = vel_bc->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  ierr = vel_bc->set_attr("valid_min", -1e6/secpera, 0); CHKERRQ(ierr); 
  ierr = vel_bc->set_attr("valid_max",  1e6/secpera, 0); CHKERRQ(ierr); 
  ierr = vel_bc->set_attr("valid_min", -1e6/secpera, 1); CHKERRQ(ierr); 
  ierr = vel_bc->set_attr("valid_max",  1e6/secpera, 1); CHKERRQ(ierr); 
  ierr = vel_bc->set_attr("_FillValue", 2e6/secpera, 0); CHKERRQ(ierr); 
  ierr = vel_bc->set_attr("_FillValue", 2e6/secpera, 1); CHKERRQ(ierr); 
  vel_bc->write_in_glaciological_units = true;
  ierr = vel_bc->set(2e6/secpera); CHKERRQ(ierr); 
  ierr = vars.add(*vel_bc); CHKERRQ(ierr);

  ierr = longitude->create(grid, "lon", true); CHKERRQ(ierr);
  ierr = longitude->set_attrs("mapping", "longitude", "degree_east", "longitude"); CHKERRQ(ierr);
  longitude->time_independent = true;
  ierr = longitude->set_attr("coordinates", ""); CHKERRQ(ierr);
  ierr = longitude->set_attr("grid_mapping", ""); CHKERRQ(ierr);
  ierr = vars.add(*longitude); CHKERRQ(ierr);

  // latitude
  ierr = latitude->create(grid, "lat", true); CHKERRQ(ierr); // has ghosts so that we can compute cell areas
  ierr = latitude->set_attrs("mapping", "latitude", "degree_north", "latitude"); CHKERRQ(ierr);
  latitude->time_independent = true;
  ierr = latitude->set_attr("coordinates", ""); CHKERRQ(ierr);
  ierr = latitude->set_attr("grid_mapping", ""); CHKERRQ(ierr);
  ierr = vars.add(*latitude); CHKERRQ(ierr);

  return 0;
}

  //! \brief Clean up.
PetscErrorCode deallocate_vars(PISMVars &variables) {
  set<string> vars = variables.keys();

  set<string>::iterator j = vars.begin();
  while (j != vars.end())
    delete variables.get(*j++);

  return 0;
}

PetscErrorCode grid_setup(IceGrid &grid) {
  PetscErrorCode ierr;
  string filename;
  PISMIO pio(&grid);
  grid_info g;

  ierr = PetscOptionsBegin(grid.com, "", "PROSS Grid options", ""); CHKERRQ(ierr);
  {
    bool flag;
    ierr = PISMOptionsString("-boot_file", "A file to bootstrap from",
                             filename, flag); CHKERRQ(ierr);
    if (!flag) {
      PetscPrintf(grid.com, "ERROR: -boot_file is required\n");
      PISMEnd();
    }

    ierr = pio.open_for_reading(filename.c_str()); CHKERRQ(ierr);
    ierr = pio.get_grid_info_2d(g); CHKERRQ(ierr);
    ierr = pio.close(); CHKERRQ(ierr); 
    grid.Mx = g.x_len;
    grid.My = g.y_len;
    grid.Mz = 2;
    grid.x0 = g.x0;
    grid.y0 = g.y0;
    grid.Lx = g.Lx;
    grid.Ly = g.Ly;

    ierr = PISMOptionsInt("-Mx", "Number of grid points in the X-direction",
                          grid.Mx, flag); CHKERRQ(ierr);
    ierr = PISMOptionsInt("-My", "Number of grid points in the Y-direction",
                          grid.My, flag); CHKERRQ(ierr);
    
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid.compute_vertical_levels(); CHKERRQ(ierr);
  ierr = grid.createDA(); CHKERRQ(ierr); 

  return 0;
}

PetscErrorCode set_surface_elevation(PISMVars &vars, const NCConfigVariable &config) {
  PetscErrorCode ierr;
  IceModelVec2S *surface, *thickness;

  surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (surface == NULL) SETERRQ(1, "surface_altitude is not available");

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  PetscReal ice_rho = config.get("ice_density"),
    ocean_rho = config.get("sea_water_density");

  ierr = thickness->copy_to(*surface); CHKERRQ(ierr);
  ierr = surface->scale(1.0 - ice_rho / ocean_rho); CHKERRQ(ierr); 

  return 0;
}


  //! \brief Read input data (by regridding).
PetscErrorCode read_input_data(IceGrid &grid, PISMVars &variables, const NCConfigVariable &config) {
  PetscErrorCode ierr;
  set<string> vars = variables.keys();
  PISMIO pio(&grid);
  grid_info g;
  string filename;

  bool boot_file_set;
  ierr = PetscOptionsBegin(grid.com, "", "PROSS I/O options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-boot_file", "A file to bootstrap from",
                             filename, boot_file_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);


  ierr = pio.open_for_reading(filename.c_str()); CHKERRQ(ierr);
  ierr = pio.get_grid_info(g); CHKERRQ(ierr);
  ierr = pio.close(); CHKERRQ(ierr);

  LocalInterpCtx lic(g, NULL, NULL, grid);

  // Read in everything except enthalpy, usurf and tauc:
  vars.erase("enthalpy");
  vars.erase("usurf");
  vars.erase("surface_altitude");
  vars.erase("tauc");
  set<string>::iterator j = vars.begin();
  while (j != vars.end()) {
    IceModelVec *var = variables.get(*j);
    ierr = var->regrid(filename.c_str(), lic, true); CHKERRQ(ierr);
    j++;
  }

  ierr = variables.get("enthalpy")->set(528668.35); CHKERRQ(ierr);
  // arbitrary; corresponds to 263.15 Kelvin at depth=0.
  // The CustomGlenIce flow law does not use it.

  // set the surface elevation:
  ierr = set_surface_elevation(variables, config); CHKERRQ(ierr); 

  // set the basal yield stress (does not matter; everything is floating)
  ierr = variables.get("tauc")->set(0.0); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode write_results(ShallowStressBalance &ssa,
                             IceGrid &grid, PISMVars &variables) {
  PetscErrorCode ierr;
  set<string> vars = variables.keys();
  string filename = "ross_new_output.nc";
  bool flag;

  ierr = PISMOptionsString("-o", "Specifies an output file",
                           filename, flag); CHKERRQ(ierr);

  PISMIO pio(&grid);

  ierr = pio.open_for_writing(filename, false, true); CHKERRQ(ierr);
  ierr = pio.append_time(0.0);
  ierr = pio.close(); CHKERRQ(ierr); 

  set<string>::iterator j = vars.begin();
  while (j != vars.end()) {
    IceModelVec *var = variables.get(*j);
    ierr = var->write(filename.c_str()); CHKERRQ(ierr);
    j++;
  }

  IceModelVec2V *vel_ssa;

  ierr = ssa.get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr);
  vel_ssa->write_in_glaciological_units = true;
  ierr = vel_ssa->write(filename.c_str()); CHKERRQ(ierr);

  IceModelVec2S cbar;

  ierr = cbar.create(grid, "cbar", false); CHKERRQ(ierr);
  ierr = cbar.set_attrs("diagnostic",
                        "magnitude of vertically-integrated horizontal velocity of ice",
                        "m s-1", ""); CHKERRQ(ierr);
  ierr = cbar.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  cbar.write_in_glaciological_units = true;

  ierr = vel_ssa->magnitude(cbar); CHKERRQ(ierr);

  ierr = cbar.write(filename.c_str()); CHKERRQ(ierr);
  
  return 0;
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);
  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  { /* This explicit scoping forces destructors to be called before PetscFinalize() */
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PROSS %s (EISMINT-Ross diagnostic velocity computation mode)\n",
                      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    string usage =
      "  pross -boot_file IN.nc -Mx number -My number [-o file.nc] [-riggs file.nc]\n"
      "where:\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -Mx         number of grid points in the x direction\n"
      "  -My         number of grid points in the y direction\n"
      "  -riggs      read RIGGS data from a file\n"
      "notes:\n"
      "  * -boot_file is required\n";

    vector<string> required;
    required.push_back("-boot_file");
    ierr = show_usage_check_req_opts(com, "pross", required, usage.c_str()); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    config.set_flag("use_ssa_velocity", true);
    config.set_flag("use_ssa_when_grounded", false);
    config.set_flag("use_constant_nuh_for_ssa", false);
    config.set("epsilon_ssa", 0.0);  // don't use this lower bound on effective viscosity

    ierr = PetscOptionsBegin(com, "", "PROSS options", ""); CHKERRQ(ierr);
    {
      bool flag;
      double rtol = config.get("ssa_relative_convergence");
      ierr = PISMOptionsReal("-ssa_rtol", "set configuration constant ssa_relative_convergence",
                             rtol, flag); CHKERRQ(ierr);
      if (flag) config.set("ssa_relative_convergence",rtol);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    IceGrid g(com, rank, size, config);

    ierr = grid_setup(g); CHKERRQ(ierr); 
    ierr = g.printInfo(1); CHKERRQ(ierr);

    PISMVars vars;

    ierr = allocate_vars(g, vars); CHKERRQ(ierr); 

    ierr = read_input_data(g, vars, config); CHKERRQ(ierr); 

    CustomGlenIce ice(g.com, "", config);

    IceBasalResistancePlasticLaw basal(config.get("plastic_regularization") / secpera, 
                                       config.get_flag("do_pseudo_plastic_till"),
                                       config.get("pseudo_plastic_q"),
                                       config.get("pseudo_plastic_uthreshold") / secpera);
    ierr = basal.printInfo(1,g.com); CHKERRQ(ierr);

    EnthalpyConverter EC(config);

    // Create the SSA solver object:
    SSAFD ssa(g, basal, ice, EC, config);

    const PetscReal
      DEFAULT_MIN_THICKNESS = 5.0, // meters
      DEFAULT_CONSTANT_HARDNESS_FOR_SSA = 1.9e8,  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
      DEFAULT_TYPICAL_STRAIN_RATE = (100.0 / secpera) / (100.0 * 1.0e3),  // typical strain rate is 100 m/yr per 
      DEFAULT_nuH = DEFAULT_MIN_THICKNESS * DEFAULT_CONSTANT_HARDNESS_FOR_SSA /
      (2.0 * pow(DEFAULT_TYPICAL_STRAIN_RATE,2./3.)); // Pa s m
    // COMPARE: 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of
    //          30 MPa yr for \bar\nu

    ssa.strength_extension->set_min_thickness(DEFAULT_MIN_THICKNESS);
    ssa.strength_extension->set_notional_strength(DEFAULT_nuH);
    ice.setHardness(DEFAULT_CONSTANT_HARDNESS_FOR_SSA);

    ierr = ssa.init(vars); CHKERRQ(ierr);

    IceModelVec2Mask *bc_mask;
    IceModelVec2V *bc_vel;
    bc_mask = dynamic_cast<IceModelVec2Mask*>(vars.get("bcflag"));
    if (bc_mask == NULL) SETERRQ(1, "bcflag is not available");

    bc_vel = dynamic_cast<IceModelVec2V*>(vars.get("velbar"));
    if (bc_vel == NULL) SETERRQ(1, "velbar is not available");

    ierr = ssa.set_boundary_conditions(*bc_mask, *bc_vel); CHKERRQ(ierr);

    ierr = verbPrintf(2,com,"* Solving the SSA stress balance ...\n"); CHKERRQ(ierr);

    ierr = ssa.update(false); CHKERRQ(ierr);

    string ssa_stdout;
    ierr = ssa.stdout_report(ssa_stdout); CHKERRQ(ierr);
    ierr = verbPrintf(3,com,ssa_stdout.c_str()); CHKERRQ(ierr);
    
    IceModelVec2V *vel_ssa;
    ierr = ssa.get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr); 

    ierr = read_riggs_and_compare(g, vars, *vel_ssa); CHKERRQ(ierr);

    ierr = compute_errors(g, vars, *vel_ssa); CHKERRQ(ierr); 

    ierr = write_results(ssa, g, vars); CHKERRQ(ierr); 
    
    ierr = deallocate_vars(vars); CHKERRQ(ierr); 
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
