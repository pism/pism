// Copyright (C) 2011 Andy Aschwanden and Constantine Khroulev
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

#include "PSElevation.hh"
#include "PISMIO.hh"


///// Elevation-dependent temperature and surface mass balance.

PetscErrorCode PSElevation::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the constant-in-time surface processes model PSElevation. Setting...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Elevation-dependent surface model options", ""); CHKERRQ(ierr);
  {
    PetscInt Nparam= 7;
    PetscReal inarray[7] = {-30, 6., 3.0, 3.0, 1500., 0., 3000.0};

    ierr = PetscOptionsGetRealArray(PETSC_NULL, "-elevation_to_artm_and_acab", inarray, &Nparam, &elev_set);
    CHKERRQ(ierr);

    T_ELA = inarray[0];
    dTdz = inarray[1]; 
    dabdz = inarray[2]; 
    dacdz = inarray[3]; 
    z_ELA = inarray[4];
    z_min = inarray[5];
    z_max = inarray[6];
    
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = verbPrintf(3, grid.com,
                    "     temperature at equilibrium line altitude T_ELA = %.2f deg C\n"
                    "     temperature gradient dTdz = %.2f (K/km)\n"
                    "     mass balance gradient in ablation area dabdz = %.2 m/year/km\n"
                    "     mass balance gradient in accumulation area dacdz = %.2f m/year/km\n"
                    "     equilibrium line altitude z_ELA = %.2f m\n"
                    "     elevation below which acab is constant z_min = %.2f m\n"
                    "     elevation above which acab is constant z_max = %.2f m\n\n",
                    T_ELA,dTdz,dabdz,dacdz,z_ELA,z_min,z_max); CHKERRQ(ierr);

  // get access to ice upper surface elevation
  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!usurf) SETERRQ(12, "ERROR: 'usurf' is not available or is wrong type in dictionary");


  // allocate IceModelVecs for storing temperature and surface mass balance fields

  acab.init_2d("acab", grid);
  acab.set_string("pism_intent", "diagnostic");
  acab.set_string("long_name",
                  "ice-equivalent surface mass balance (accumulation/ablation) rate");
  acab.set_string("standard_name",
                  "land_ice_surface_specific_mass_balance");
  ierr = acab.set_units("m s-1"); CHKERRQ(ierr);
  ierr = acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  artm.init_2d("artm", grid);
  artm.set_string("pism_intent", "diagnostic");
  artm.set_string("long_name",
                  "ice temperature at the ice surface");
  ierr = artm.set_units("K"); CHKERRQ(ierr);

  // parameterizing the ice surface temperature 'artm'
  ierr = verbPrintf(2, grid.com,"    parameterizing the ice surface temperature 'artm' ... \n"); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com, 
                    "      ice temperature at the ice surface (artm) is piecewise-linear function of surface altitude (usurf):\n"
                    "          artm = %5.2f K + %.2f K/km * (usurf - z_ELA)\n",
                    T_ELA+273.15 , dTdz); CHKERRQ(ierr);

  // parameterizing the ice surface mass balance 'acab'
  ierr = verbPrintf(2, grid.com,"    parameterizing the ice surface mass balance 'acab' ... \n"); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
                    "      surface mass balance (acab) is piecewise-linear function of surface altitue (usurf):\n"
                    "                /  %5.2f m/a                    for            usurf < %3.f m\n"
                    "          acab =  %5.2f m/a/km * (usurf-%3.f m) for   %3.f m < usurf < %3.f m\n"
                    "               \\  %5.2f m/a/km * (usurf-%3.f m) for   %3.f m < usurf < %3.f m\n"
                    "                \\ %5.2f m/a                     for   %3.f m < usurf\n",
                    dabdz / 1000 * (z_min - z_ELA), z_min,
                    dabdz, z_ELA, z_min, z_ELA, 
                    dacdz, z_ELA, z_ELA, z_max,
                    dacdz / 1000 * (z_max - z_ELA), z_max); CHKERRQ(ierr);


  return 0;
}

PetscErrorCode PSElevation::ice_surface_mass_flux(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						 IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal z_minELA = z_min - z_ELA;
  PetscReal z_maxELA = z_max - z_ELA;
  string history  = "elevation-dependent surface mass balance\n";

  ierr = result.begin_access();   CHKERRQ(ierr);
  ierr = usurf->begin_access();   CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {      
      PetscReal z = (*usurf)(i,j);
      if (z < z_min)
        {
          result(i,j) = dabdz/1000/secpera * z_minELA;
        }
      else if ((z >= z_min) && (z < z_ELA))
        {
          result(i,j) = dabdz/1000/secpera * (z - z_ELA);
        }
      else if ((z >= z_ELA) && (z <= z_max))
        {
          result(i,j) = dacdz/1000/secpera * (z - z_ELA);
        }
      else if (z > z_max)
        {
          result(i,j) = dacdz/1000/secpera * z_maxELA;
        }
      else
        {
          SETERRQ(1,"HOW DID I GET HERE? ... ending...\n");
        }
      ierr = verbPrintf(5, grid.com,"!!!!! z=%f, acab=%f\n",z,result(i,j)); CHKERRQ(ierr);
    }
  }
  ierr = usurf->end_access();   CHKERRQ(ierr);
  ierr = result.end_access();   CHKERRQ(ierr);

  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSElevation::ice_surface_temperature(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						   IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "elevation-dependent ice surface temperature \n";

  ierr = result.begin_access();   CHKERRQ(ierr);
  ierr = usurf->begin_access();   CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {      
      PetscReal z = (*usurf)(i,j);
      result(i,j) = 273.15 + T_ELA + dTdz /1000 * (z - z_ELA);
      ierr = verbPrintf(5, grid.com,"!!!!! z=%f, artm=%f\n",z,result(i,j)); CHKERRQ(ierr);
    }
  }
  ierr = usurf->end_access();   CHKERRQ(ierr);
  ierr = result.end_access();   CHKERRQ(ierr);

  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

void PSElevation::add_vars_to_output(string keyword, set<string> &result) {
  if (keyword != "small") {
    result.insert("artm");
    result.insert("acab");
  }
}

PetscErrorCode PSElevation::define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;
  int varid;

  ierr = PISMSurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  if (set_contains(vars, "artm")) {
    ierr = artm.define(nc, varid, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    ierr = acab.define(nc, varid, nctype, true); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSElevation::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "artm", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(artm, 0); CHKERRQ(ierr);

    ierr = ice_surface_temperature(t, dt, tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "acab", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(acab, 0); CHKERRQ(ierr);

    ierr = ice_surface_mass_flux(t, dt, tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}
