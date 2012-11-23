// Copyright (C) 2010, 2011, 2012 Ed Bueler and Constantine Khroulev
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

static char help[] =
  "\nSIAFD_TEST\n"
  "  Testing program for SIA, time-independent calculations separate from\n"
  "  IceModel. Uses verification test F. Also may be used in a PISM software"
  "(regression) test.\n\n";

#include "pism_const.hh"
#include "pism_options.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "PIO.hh"
#include "NCVariable.hh"
#include "PISMStressBalance.hh"
#include "SIAFD.hh"
#include "exactTestsFG.h"
#include "basal_resistance.hh"
#include "enthalpyConverter.hh"
#include "SSB_Modifier.hh"
#include "ShallowStressBalance.hh"
#include "PISMVars.hh"

PetscErrorCode computeSigmaErrors(const NCConfigVariable &config,
                                  IceModelVec3 &Sigma,
                                  IceModelVec2S &thickness,
                                  IceGrid &grid,
                                  PetscScalar &gmaxSigmaerr,
                                  PetscScalar &gavSigmaerr) {

  PetscErrorCode ierr;
  PetscScalar    maxSigerr = 0.0, avSigerr = 0.0, avcount = 0.0;
  PetscScalar    **H;
  const PetscInt Mz = grid.Mz;
  
  const PetscScalar LforFG = 750000; // m

  const PetscScalar
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity");

  PetscScalar   *dummy1, *dummy2, *dummy3, *dummy4, *Sigex;
  PetscScalar   junk0, junk1;
  
  Sigex = new PetscScalar[Mz];
  dummy1 = new PetscScalar[Mz];  dummy2 = new PetscScalar[Mz];
  dummy3 = new PetscScalar[Mz];  dummy4 = new PetscScalar[Mz];

  PetscScalar *Sig;

  ierr = thickness.get_array(H); CHKERRQ(ierr);
  ierr = Sigma.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar xx = grid.x[i], yy = grid.y[j],
        r = sqrt(PetscSqr(xx) + PetscSqr(yy));
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        bothexact(0.0,r,&grid.zlevels[0],Mz,0.0,
                  &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,dummy4);

        for (PetscInt k = 0; k < Mz; k++)
          Sigex[k] = Sigex[k] * ice_rho * ice_c; // scale exact Sigma to J/(s m^3)
        const PetscInt ks = grid.kBelowHeight(H[i][j]);
        ierr = Sigma.getInternalColumn(i,j,&Sig); CHKERRQ(ierr);
        for (PetscInt k = 0; k < ks; k++) {  // only eval error if below num surface
          const PetscScalar Sigerr = PetscAbs(Sig[k] - Sigex[k]);
          maxSigerr = PetscMax(maxSigerr,Sigerr);
          avcount += 1.0;
          avSigerr += Sigerr;
        }
      }
    }
  }
  ierr = thickness.end_access(); CHKERRQ(ierr);
  ierr = Sigma.end_access(); CHKERRQ(ierr);

  delete [] Sigex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  
  ierr = PISMGlobalMax(&maxSigerr, &gmaxSigmaerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avSigerr, &gavSigmaerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PISMGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavSigmaerr = gavSigmaerr/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode computeSurfaceVelocityErrors(IceGrid &grid,
                                            IceModelVec2S &vH,
                                            IceModelVec3 &u3,
                                            IceModelVec3 &v3,
                                            IceModelVec3 &w3,
                                            PetscScalar &gmaxUerr,
                                            PetscScalar &gavUerr,
                                            PetscScalar &gmaxWerr,
                                            PetscScalar &gavWerr) {

  PetscErrorCode ierr;
  PetscScalar    maxUerr = 0.0, maxWerr = 0.0, avUerr = 0.0, avWerr = 0.0;
  PetscScalar    **H;

  const PetscScalar LforFG = 750000; // m

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar xx = grid.x[i], yy = grid.y[j],
        r = sqrt(PetscSqr(xx) + PetscSqr(yy));
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        PetscScalar radialUex,wex;
        PetscScalar dummy0,dummy1,dummy2,dummy3,dummy4;
        bothexact(0.0,r,&(H[i][j]),1,0.0,
                  &dummy0,&dummy1,&dummy2,&radialUex,&wex,&dummy3,&dummy4);

        const PetscScalar uex = (xx/r) * radialUex;
        const PetscScalar vex = (yy/r) * radialUex;
        // note that because getValZ does linear interpolation and H[i][j] is not exactly at
        // a grid point, this causes nonzero errors even with option -eo
        const PetscScalar Uerr = sqrt(PetscSqr(u3.getValZ(i,j,H[i][j]) - uex)
                                      + PetscSqr(v3.getValZ(i,j,H[i][j]) - vex));
        maxUerr = PetscMax(maxUerr,Uerr);
        avUerr += Uerr;
        const PetscScalar Werr = PetscAbs(w3.getValZ(i,j,H[i][j]) - wex);
        maxWerr = PetscMax(maxWerr,Werr);
        avWerr += Werr;
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&maxUerr, &gmaxUerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&maxWerr, &gmaxWerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avUerr, &gavUerr, grid.com); CHKERRQ(ierr);
  gavUerr = gavUerr/(grid.Mx*grid.My);
  ierr = PISMGlobalSum(&avWerr, &gavWerr, grid.com); CHKERRQ(ierr);
  gavWerr = gavWerr/(grid.Mx*grid.My);
  return 0;
}


PetscErrorCode enthalpy_from_temperature_cold(EnthalpyConverter &EC,
                                              IceGrid &grid,
                                              IceModelVec2S &thickness,
                                              IceModelVec3 &temperature,
                                              IceModelVec3 &enthalpy) {
  PetscErrorCode ierr;
  
  ierr = temperature.begin_access(); CHKERRQ(ierr);
  ierr = enthalpy.begin_access(); CHKERRQ(ierr);
  ierr = thickness.begin_access(); CHKERRQ(ierr);

  PetscScalar *T_ij, *E_ij; // columns of these values
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = temperature.getInternalColumn(i,j,&T_ij); CHKERRQ(ierr);
      ierr = enthalpy.getInternalColumn(i,j,&E_ij); CHKERRQ(ierr);

      for (PetscInt k=0; k<grid.Mz; ++k) {
        PetscReal depth = thickness(i,j) - grid.zlevels[k];
        ierr = EC.getEnthPermissive(T_ij[k], 0.0,
                                    EC.getPressureFromDepth(depth),
                                    E_ij[k]); CHKERRQ(ierr);
      }

    }
  }

  ierr = enthalpy.end_access(); CHKERRQ(ierr);
  ierr = temperature.end_access(); CHKERRQ(ierr);
  ierr = thickness.end_access(); CHKERRQ(ierr);

  ierr = enthalpy.beginGhostComm(); CHKERRQ(ierr);
  ierr = enthalpy.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


//! \brief Set the test F initial state.
PetscErrorCode setInitStateF(IceGrid &grid,
                             EnthalpyConverter &EC,
                             IceModelVec2S *bed,
                             IceModelVec2Int *mask,
                             IceModelVec2S *surface,
                             IceModelVec2S *thickness,
                             IceModelVec3 *enthalpy) {
  PetscErrorCode ierr;
  PetscInt        Mz=grid.Mz;
  PetscScalar     **H;
  PetscScalar     *dummy1, *dummy2, *dummy3, *dummy4, *dummy5;

  PetscReal ST = 1.67e-5,
    Tmin = 223.15,  // K
    LforFG = 750000; // m

  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];
  dummy5=new PetscScalar[Mz];

  ierr = bed->set(0); CHKERRQ(ierr);
  ierr = mask->set(MASK_GROUNDED); CHKERRQ(ierr);

  PetscScalar *T;
  T = new PetscScalar[grid.Mz];

  ierr = thickness->get_array(H); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r = grid.radius(i,j),
        Ts = Tmin + ST * r;
      if (r > LforFG - 1.0) { // if (essentially) outside of sheet
        H[i][j] = 0.0;
        for (PetscInt k = 0; k < Mz; k++)
          T[k]=Ts;
      } else {
        r = PetscMax(r,1.0); // avoid singularity at origin
        bothexact(0.0,r,&grid.zlevels[0],Mz,0.0,
                  &H[i][j],dummy5,T,dummy1,dummy2,dummy3,dummy4);
      }
      ierr = enthalpy->setInternalColumn(i,j,T); CHKERRQ(ierr);
    }
  }

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr =  enthalpy->end_access(); CHKERRQ(ierr);

  ierr = thickness->beginGhostComm(); CHKERRQ(ierr);
  ierr = thickness->endGhostComm(); CHKERRQ(ierr);

  ierr = thickness->copy_to(*surface); CHKERRQ(ierr);

  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  delete [] T; delete [] dummy5;

  ierr = enthalpy_from_temperature_cold(EC, grid, *thickness,
                                        *enthalpy,
                                        *enthalpy); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode reportErrors(const NCConfigVariable &config,
                            IceGrid &grid,
                            IceModelVec2S *thickness,
                            IceModelVec3 *u_sia, IceModelVec3 *v_sia,
                            IceModelVec3 *w_sia,
                            IceModelVec3 *Sigma
                            ) {
  PetscErrorCode ierr;

  // Sigma errors if appropriate; reported in 10^6 J/(s m^3)
  PetscScalar maxSigerr, avSigerr;
  ierr = computeSigmaErrors(config, *Sigma, *thickness,
                            grid,
                            maxSigerr, avSigerr); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com, 
                    "Sigma     :      maxSig       avSig\n"); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f\n", 
                    maxSigerr*1.0e6, avSigerr*1.0e6); CHKERRQ(ierr);

  // surface velocity errors if exact values are available; reported in m/a
  PetscScalar maxUerr, avUerr, maxWerr, avWerr;
  ierr = computeSurfaceVelocityErrors(grid, *thickness,
                                      *u_sia,
                                      *v_sia,
                                      *w_sia,
                                      maxUerr, avUerr,
                                      maxWerr, avWerr); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com, 
                    "surf vels :     maxUvec      avUvec        maxW         avW\n"); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n", 
                    convert(maxUerr, "m/s", "m/year"),
                    convert(avUerr, "m/s", "m/year"),
                    convert(maxWerr, "m/s", "m/year"),
                    convert(avWerr, "m/s", "m/year")); CHKERRQ(ierr);

  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
  
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {  
    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    config.set_flag("compute_grain_size_using_age", false);

    PetscBool usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage of SIAFD_TEST:\n"
                  "  run siafd_test -Mx <number> -My <number> -Mz <number> -o foo.nc\n"
                  "\n");
    }

    IceGrid grid(com, rank, size, config);

    grid.Lx = grid.Ly = 900e3;
    grid.Lz = 4000;
    grid.Mx = grid.My = 61;
    grid.Mz = 61;
    
    string output_file = "siafd_test_F.nc";
    ierr = PetscOptionsBegin(grid.com, "", "SIAFD_TEST options", ""); CHKERRQ(ierr);
    {
      bool flag;
      ierr = PISMOptionsInt("-Mx", "Number of grid points in the X direction",
                            grid.Mx, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-My", "Number of grid points in the X direction",
                            grid.My, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-Mz", "Number of vertical grid levels",
                            grid.Mz, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Set the output file name",
                               output_file, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_vertical_levels(); CHKERRQ(ierr); 
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);
    ierr = grid.printInfo(1); CHKERRQ(ierr);
    //ierr = grid.printVertLevels(1); CHKERRQ(ierr); 

    ICMEnthalpyConverter EC(config);
    ThermoGlenArrIce ice(grid.com, "", config, &EC);

    IceModelVec2S vh, vH, vbed;
    IceModelVec2Int vMask;
    IceModelVec3 enthalpy,
      age;                      // is not used (and need not be allocated)
    const PetscInt WIDE_STENCIL = grid.max_stencil_width;

    PISMVars vars;

    // ice upper surface elevation
    ierr = vh.create(grid, "usurf", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vh.set_attrs("diagnostic", "ice upper surface elevation",
          "m", "surface_altitude"); CHKERRQ(ierr);
    ierr = vars.add(vh); CHKERRQ(ierr);

    // land ice thickness
    ierr = vH.create(grid, "thk", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vH.set_attrs("model_state", "land ice thickness",
          "m", "land_ice_thickness"); CHKERRQ(ierr);
    ierr = vH.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = vars.add(vH); CHKERRQ(ierr);

    // bedrock surface elevation
    ierr = vbed.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vbed.set_attrs("model_state", "bedrock surface elevation",
          "m", "bedrock_altitude"); CHKERRQ(ierr);
    ierr = vars.add(vbed); CHKERRQ(ierr);

    // age of the ice; is not used here
    ierr = age.create(grid, "age", false); CHKERRQ(ierr);
    ierr = age.set_attrs("diagnostic", "age of the ice", "s", ""); CHKERRQ(ierr);
    ierr = age.set_glaciological_units("year"); CHKERRQ(ierr);
    age.write_in_glaciological_units = true;
    ierr = vars.add(age); CHKERRQ(ierr);

    // enthalpy in the ice
    ierr = enthalpy.create(grid, "enthalpy", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = enthalpy.set_attrs("model_state",
                              "ice enthalpy (includes sensible heat, latent heat, pressure)",
                              "J kg-1", ""); CHKERRQ(ierr);
    ierr = vars.add(enthalpy); CHKERRQ(ierr);

    // grounded_dragging_floating integer mask
    ierr = vMask.create(grid, "mask", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vMask.set_attrs("model_state", "grounded_dragging_floating integer mask",
			 "", ""); CHKERRQ(ierr);
    vector<double> mask_values(4);
    mask_values[0] = MASK_ICE_FREE_BEDROCK;
    mask_values[1] = MASK_GROUNDED;
    mask_values[2] = MASK_FLOATING;
    mask_values[3] = MASK_ICE_FREE_OCEAN;
    ierr = vMask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
    ierr = vMask.set_attr("flag_meanings",
			"ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
		  CHKERRQ(ierr);
    vMask.output_data_type = PISM_BYTE;
    ierr = vars.add(vMask); CHKERRQ(ierr);

    // This is never used (but it is a required argument of the
    // PISMStressBalance constructor).
    IceBasalResistancePlasticLaw basal(
           config.get("plastic_regularization", "1/year", "1/second"), 
           config.get_flag("do_pseudo_plastic_till"),
           config.get("pseudo_plastic_q"),
           config.get("pseudo_plastic_uthreshold", "m/year", "m/second"));

    // Create the SIA solver object:

    // We use SIA_Nonsliding and not SIAFD here because we need the z-component
    // of the ice velocity, which is computed using incompressibility of ice in
    // PISMStressBalance::compute_vertical_velocity().
    SIAFD *sia = new SIAFD(grid, EC, config);
    SSB_Trivial *trivial_stress_balance = new SSB_Trivial(grid, basal, EC, config);

    PISMStressBalance stress_balance(grid,
                                     trivial_stress_balance, sia,
                                     NULL, // no ocean model
                                     config);

    // fill the fields:
    ierr = setInitStateF(grid, EC,
                         &vbed, &vMask, &vh, &vH,
                         &enthalpy); CHKERRQ(ierr);

    // Allocate the SIA solver:
    ierr = stress_balance.init(vars); CHKERRQ(ierr);

    // Solve (fast==true means "no 3D update and no strain heating computation"):
    bool fast = false;
    ierr = stress_balance.update(fast); CHKERRQ(ierr); 

    // Report errors relative to the exact solution:
    IceModelVec3 *u_sia, *v_sia, *w_sia, *sigma;
    ierr = stress_balance.get_3d_velocity(u_sia, v_sia, w_sia); CHKERRQ(ierr); 

    ierr = stress_balance.get_volumetric_strain_heating(sigma); CHKERRQ(ierr); 

    ierr = reportErrors(config, grid,
                        &vH, u_sia, v_sia, w_sia, sigma); CHKERRQ(ierr);

    // Write results to an output file:
    PIO pio(grid, "guess_format");

    ierr = pio.open(output_file, PISM_WRITE); CHKERRQ(ierr);
    ierr = pio.def_time(config.get_string("time_dimension_name"),
                        config.get_string("calendar"),
                        grid.time->CF_units()); CHKERRQ(ierr);
    ierr = pio.append_time(config.get_string("time_dimension_name"), 0.0);
    ierr = pio.close(); CHKERRQ(ierr);

    ierr = vh.write(output_file); CHKERRQ(ierr);
    ierr = vH.write(output_file); CHKERRQ(ierr);
    ierr = vMask.write(output_file); CHKERRQ(ierr);
    ierr = vbed.write(output_file); CHKERRQ(ierr);
    // ierr = enthalpy.write(output_file); CHKERRQ(ierr);
    ierr = u_sia->write(output_file); CHKERRQ(ierr);
    ierr = v_sia->write(output_file); CHKERRQ(ierr);
    ierr = w_sia->write(output_file); CHKERRQ(ierr);
    ierr = sigma->write(output_file); CHKERRQ(ierr);
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
