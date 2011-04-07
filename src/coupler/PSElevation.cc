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
    PetscInt Nartmparam= 4;
    PetscReal artmarray[4] = {-5,0,1325,1350};

    ierr = PetscOptionsGetRealArray(PETSC_NULL, "-artm", artmarray, &Nartmparam, &elev_artm_set);
    CHKERRQ(ierr);

    artm_min = artmarray[0];
    artm_max = artmarray[1];
    z_artm_min = artmarray[2];
    z_artm_max = artmarray[3];

    PetscInt Nacabparam= 5;
    PetscReal acabarray[5] = {-3,3,1100,1500,1700};

    ierr = PetscOptionsGetRealArray(PETSC_NULL, "-acab", acabarray, &Nacabparam, &elev_acab_set);
    CHKERRQ(ierr);

    acab_min = acabarray[0];
    acab_max = acabarray[1];
    z_acab_min = acabarray[2];
    z_ELA = acabarray[3];
    z_acab_max = acabarray[4];

    PetscInt Nlimitsparam= 2;
    PetscReal limitsarray[2] = {0,0};

    ierr = PetscOptionsGetRealArray(PETSC_NULL, "-acab_limits", limitsarray, &Nlimitsparam, &acab_limits_set);
    CHKERRQ(ierr);

    if (acab_limits_set)
      {
        acab_limit_min = limitsarray[0];
        acab_limit_max = limitsarray[1];
      }
    else
      {
        acab_limit_min = acab_min;
        acab_limit_max = acab_max;
      }

  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = verbPrintf(3, grid.com,
                    "     temperature at %.0f m a.s.l. = %.2f deg C\n"
                    "     temperature at %.0f m a.s.l. = %.2f deg C\n"
                    "     mass balance at and below %.0f m a.s.l. = %.2f m/a\n"
                    "     mass balance at and below %.0f m a.s.l. = %.2f m/a\n"
                    "     equilibrium line altitude z_ELA = %.2f m a.s.l.\n",
                    z_artm_min, artm_min, z_artm_max, artm_max, z_acab_min, acab_min, z_acab_max, acab_max, z_ELA); CHKERRQ(ierr);

  // get access to ice upper surface elevation
  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!usurf) SETERRQ(12, "ERROR: 'usurf' is not available or is wrong type in dictionary");


  // allocate NCSpatialVariables for storing temperature and surface mass balance fields

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
                    "                 /  %2.2f K                            for            usurf < %.f m\n"
                    "         artm = |   %5.2f K + %5.3f * (usurf - %.f m) for   %.f m < usurf < %.f m\n"
                    "                 \\  %5.2f K                            for   %.f m < usurf\n",
                    artm_min + 273.15, z_artm_min,
                    artm_min + 273.15, (artm_max-artm_min)/(z_artm_max-z_artm_min), z_artm_min, z_artm_min, z_artm_max,
                    artm_max + 273.13, z_artm_max); CHKERRQ(ierr);

  // parameterizing the ice surface mass balance 'acab'
  ierr = verbPrintf(2, grid.com,"    parameterizing the ice surface mass balance 'acab' ... \n"); CHKERRQ(ierr);

  if (acab_limits_set)
    {
      ierr = verbPrintf(2, grid.com,"    option '-acab_limits' seen, limiting upper and lower bounds ... \n"); CHKERRQ(ierr);
    }

  ierr = verbPrintf(2, grid.com, 
                    "      surface mass balance (acab) is piecewise-linear function of surface altitue (usurf):\n"
                    "                  /  %5.2f m/a                       for          usurf < %3.f m\n"
                    "          acab = |    %5.3f 1/a * (usurf-%.0f m)     for %3.f m < usurf < %3.f m\n"
                    "                  \\   %5.3f 1/a * (usurf-%.0f m)     for %3.f m < usurf < %3.f m\n"
                    "                   \\ %5.2f m/a                       for %3.f m < usurf\n",
                    acab_limit_min, z_acab_min,
                    -acab_min/(z_ELA - z_acab_min),z_acab_min, z_acab_min, z_ELA, 
                    acab_max/(z_acab_max - z_ELA),z_acab_max, z_ELA, z_acab_max, 
                    acab_limit_max, z_acab_max); CHKERRQ(ierr);


  return 0;
}

PetscErrorCode PSElevation::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal dabdz = -acab_min/(z_ELA - z_acab_min);
  PetscReal dacdz = acab_max/(z_acab_max - z_ELA);
  string history  = "elevation-dependent surface mass balance\n";

  ierr = result.begin_access();   CHKERRQ(ierr);
  ierr = usurf->begin_access();   CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {      
      PetscReal z = (*usurf)(i,j);
      if (z < z_acab_min)
        {
          result(i,j) = acab_limit_min / secpera;
        }
      else if ((z >= z_acab_min) && (z < z_ELA))
        {
          result(i,j) = dabdz / secpera * (z - z_ELA);
        }
      else if ((z >= z_ELA) && (z <= z_acab_max))
        {
          result(i,j) = dacdz / secpera * (z - z_ELA);
        }
      else if (z > z_acab_max)
        {
          result(i,j) = acab_limit_max / secpera;
        }
      else
        {
          SETERRQ(1,"HOW DID I GET HERE? ... ending...\n");
        }
      ierr = verbPrintf(5, grid.com,"!!!!! z=%.2f, acab=%.2f\n",z,result(i,j)*secpera); CHKERRQ(ierr);
    }
  }
  ierr = usurf->end_access();   CHKERRQ(ierr);
  ierr = result.end_access();   CHKERRQ(ierr);

  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSElevation::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "elevation-dependent ice surface temperature \n";

  ierr = result.begin_access();   CHKERRQ(ierr);
  ierr = usurf->begin_access();   CHKERRQ(ierr);
  PetscReal dTdz = (artm_max - artm_min)/(z_artm_max - z_artm_min);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {      
      PetscReal z = (*usurf)(i,j);
      if (z <= z_artm_min)
        {
          result(i,j) = 273.15 + artm_min;
        }
      else if ((z > z_artm_min) && (z < z_artm_max))
        {
          result(i,j) = 273.15 + artm_min + dTdz * (z - z_artm_min);    
        }
      else if (z >= z_artm_max)
        {
          result(i,j) = 273.15 + artm_max;
        }
      else
        {
          SETERRQ(1,"HOW DID I GET HERE? ... ending...\n");
        }        
      ierr = verbPrintf(5, grid.com,"!!!!! z=%.2f, artm_min=%.2f, dTdz=%.2f, dz=%.2f, artm=%.2f\n",z,artm_min+273.15,dTdz,z - z_artm_min,result(i,j)); CHKERRQ(ierr);
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

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "acab", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(acab, 0); CHKERRQ(ierr);

    ierr = ice_surface_mass_flux(tmp); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}
