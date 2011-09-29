// Copyright (C) 2011 Ed Bueler
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
  "Variable k(T) version of pismr.  Conductivity depends additionally on\n"
  "temperature (enthalpy) in cold ice.\n";

#include <petsc.h>
#include "IceGrid.hh"
#include "iceModel.hh"

#include "varkenthSystem.hh"
#include "DrainageCalculator.hh"
#include "Mask.hh"

class IcevarkModel : public IceModel {
public:
  IcevarkModel(IceGrid &g, NCConfigVariable &c, NCConfigVariable &o)
     : IceModel(g,c,o) {};
protected:
  virtual PetscErrorCode enthalpyAndDrainageStep(
                      PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                      PetscScalar* bulgeCount);
};


PetscErrorCode IcevarkModel::enthalpyAndDrainageStep(
                      PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                      PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  ierr = verbPrintf(3,grid.com,
        "\n IcevarkModel::enthalpyAndDrainageStep() called \n"); CHKERRQ(ierr);

  if (config.get_flag("do_cold_ice_methods")) {
    SETERRQ(1,
      "PISM ERROR:  IcevarkModel::enthalpyAndDrainageStep() called but do_cold_ice_methods==true\n");
  }

  const PetscReal dt_secs = dt_years_TempAge * secpera;

  // get fine grid levels in ice
  PetscInt    fMz = grid.Mz_fine;  
  PetscScalar fdz = grid.dz_fine;
  vector<double> &fzlev = grid.zlevels_fine;

  const PetscScalar
    p_air     = config.get("surface_pressure"),
    ice_k     = config.get("ice_thermal_conductivity"),
    ice_c     = config.get("ice_specific_heat_capacity"),
    ice_K     = ice_k / ice_c, // enthalpy-conductivity for cold ice
    L         = config.get("water_latent_heat_fusion"),  // J kg-1
    bulgeEnthMax  = config.get("enthalpy_cold_bulge_max"), // J kg-1
    bwat_decay_rate = config.get("bwat_decay_rate"),   // m s-1
    bwat_max = config.get("bwat_max");                 // m

  DrainageCalculator dc(config);
  
  IceModelVec2S *Rb;
  IceModelVec3 *u3, *v3, *w3, *Sigma3;
  ierr = stress_balance->get_basal_frictional_heating(Rb); CHKERRQ(ierr);
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);
  ierr = stress_balance->get_volumetric_strain_heating(Sigma3); CHKERRQ(ierr); 

  PetscScalar *Enthnew;
  Enthnew = new PetscScalar[fMz];  // new enthalpy in column

  varkenthSystemCtx esys(config, Enth3, fMz, "varkenth");
  ierr = esys.initAllColumns(grid.dx, grid.dy, dt_secs, fdz); CHKERRQ(ierr);

  bool viewOneColumn;
  ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

  if (getVerbosityLevel() >= 4) {  // view: all column-independent constants correct?
    ierr = EC->viewConstants(NULL); CHKERRQ(ierr);
    ierr = esys.viewConstants(NULL, false); CHKERRQ(ierr);
  }

  // now get map-plane coupler fields: Dirichlet upper surface boundary and
  //    mass balance lower boundary under shelves
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(artm);
    ierr = surface->ice_surface_liquid_water_fraction(liqfrac_surface); CHKERRQ(ierr);
    CHKERRQ(ierr);
  } else {
    SETERRQ(4,"PISM ERROR: surface == PETSC_NULL");
  }
  if (ocean != PETSC_NULL) {
    ierr = ocean->shelf_base_mass_flux(shelfbmassflux);
        CHKERRQ(ierr);
    ierr = ocean->shelf_base_temperature(shelfbtemp);
        CHKERRQ(ierr);
  } else {
    SETERRQ(5,"PISM ERROR: ocean == PETSC_NULL");
  }

  IceModelVec2S G0 = vWork2d[0];
  ierr = G0.set_attrs("internal","upward geothermal flux at z=0","W m-2", ""); CHKERRQ(ierr);
  ierr = G0.set_glaciological_units("mW m-2");
  if (btu) {
    ierr = btu->get_upward_geothermal_flux(G0); CHKERRQ(ierr);
  } else {
    SETERRQ(3,"PISM ERROR: PISMBedThermalUnit* btu == PETSC_NULL in enthalpyAndDrainageStep()");
  }

  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

  // get other map-plane fields
  ierr = liqfrac_surface.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vbwat.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  ierr = Rb->begin_access(); CHKERRQ(ierr);
  ierr = G0.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  // these are accessed a column at a time
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = Sigma3->begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vWork3d.begin_access(); CHKERRQ(ierr);

  PetscInt liquifiedCount = 0;

  MaskQuery mask(vMask);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // for fine grid; this should *not* be replaced by call to grid.kBelowHeight()
      const PetscInt ks = static_cast<PetscInt>(floor(vH(i,j)/fdz));
#ifdef PISM_DEBUG
      // check if ks is valid
      if ((ks < 0) || (ks >= grid.Mz_fine)) {
        PetscPrintf(grid.com,
                    "ERROR: ks = %d computed at i = %d, j = %d is invalid,"
                    " possibly because of invalid ice thickness.\n",
                    ks, i, j);
        SETERRQ(1, "invalid ks");
      }
#endif

      const bool ice_free_column = (ks == 0),
                 is_floating     = mask.ocean(i,j);

      // enthalpy and pressures at top of ice
      const PetscScalar p_ks = EC->getPressureFromDepth(vH(i,j) - fzlev[ks]); // FIXME task #7297
      PetscScalar Enth_ks;
      ierr = EC->getEnthPermissive(artm(i,j), liqfrac_surface(i,j), p_ks, Enth_ks); CHKERRQ(ierr);

      // deal completely with columns with no ice; enthalpy, vbwat, vbmr all need setting
      if (ice_free_column) {
        ierr = vWork3d.setColumn(i,j,Enth_ks); CHKERRQ(ierr);
        if (mask.floating_ice(i,j)) {
          // if floating then assume-maximally saturated till to avoid "shock"
          //   when grounding line advances
          vbwat(i,j) = bwat_max;
          vbmr(i,j) = shelfbmassflux(i,j);
        } else {
          // either truely no ice or grounded or both; either way zero-out subglacial fields
          vbwat(i,j) = 0.0;  // no stored water on ice free land
          vbmr(i,j) = 0.0;    // no basal melt rate; melting is a surface process
                              //   on ice free land
        }

        goto donewithcolumn;
      } // end of if (ice_free_column)

      { // explicit scoping to deal with goto and initializers

        // ignore advection and strain heating in ice if isMarginal
        const bool isMarginal = checkThinNeigh(
                                 vH(i+1,j),vH(i+1,j+1),vH(i,j+1),vH(i-1,j+1),
                                 vH(i-1,j),vH(i-1,j-1),vH(i,j-1),vH(i+1,j-1)  );

        ierr = Enth3.getValColumn(i,j,ks,esys.Enth); CHKERRQ(ierr);
        ierr = w3->getValColumn(i,j,ks,esys.w); CHKERRQ(ierr);

        PetscScalar lambda;
        ierr = getEnthalpyCTSColumn(p_air, vH(i,j), ks, esys.Enth, esys.w, // FIXME task #7297
                                    &lambda, &esys.Enth_s); CHKERRQ(ierr);
        if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1

        // if there is subglacial water, don't allow ice base enthalpy to be below
        // pressure-melting; that is, assume subglacial water is at the pressure-
        // melting temperature and enforce continuity of temperature
        if ((vbwat(i,j) > 0.0) && (esys.Enth[0] < esys.Enth_s[0])) { 
          esys.Enth[0] = esys.Enth_s[0];
        }

        const bool base_is_cold = (esys.Enth[0] < esys.Enth_s[0]);
        const PetscScalar p1 = EC->getPressureFromDepth(vH(i,j) - fdz); // FIXME task #7297
        const bool k1_istemperate = EC->isTemperate(esys.Enth[1], p1); // level  z = + \Delta z

        // can now determine melt, but only preliminarily because of drainage,
        //   from heat flux out of bedrock, heat flux into ice, and frictional heating
        if (is_floating) {
          vbmr(i,j) = shelfbmassflux(i,j);
        } else {
          if (base_is_cold) {
              vbmr(i,j) = 0.0;  // zero melt rate if cold base
          } else {
            PetscScalar hf_up;
            if (k1_istemperate) {
              const PetscScalar pbasal = EC->getPressureFromDepth(vH(i,j)); // FIXME task #7297
              hf_up = - ice->k * (EC->getMeltingTemp(p1) - EC->getMeltingTemp(pbasal)) / fdz;
            } else {
              hf_up = - ice_K * (esys.Enth[1] - esys.Enth[0]) / fdz;
            }
            // compute basal melt rate from flux balance; vbmr = - Mb / rho in
            //   efgis paper; after we compute it we make sure there is no
            //   refreeze if there is no available basal water
            vbmr(i,j) = ( (*Rb)(i,j) + G0(i,j) - hf_up ) / (ice->rho * L);
            if ((vbwat(i,j) <= 0) && (vbmr(i,j) < 0))    vbmr(i,j) = 0.0;
          }
        }

        // now set-up for solve in ice; note esys.Enth[], esys.w[],
        //   esys.Enth_s[] are already filled
        ierr = esys.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

        ierr = u3->getValColumn(i,j,ks,esys.u); CHKERRQ(ierr);
        ierr = v3->getValColumn(i,j,ks,esys.v); CHKERRQ(ierr);
        ierr = Sigma3->getValColumn(i,j,ks,esys.Sigma); CHKERRQ(ierr);

        ierr = esys.setSchemeParamsThisColumn(isMarginal, lambda); CHKERRQ(ierr);
        ierr = esys.setBoundaryValuesThisColumn(Enth_ks); CHKERRQ(ierr);

        // determine lowest-level equation at bottom of ice; see decision chart
        //   in [\ref AschwandenBuelerKhroulevBlatter], and page documenting BOMBPROOF
        if (is_floating) {
          // floating base: Dirichlet application of known temperature from ocean
          //   coupler; assumes base of ice shelf has zero liquid fraction
          PetscScalar Enth0;
          ierr = EC->getEnthPermissive(shelfbtemp(i,j), 0.0, EC->getPressureFromDepth(vH(i,j)),
                                       Enth0); CHKERRQ(ierr);
          ierr = esys.setDirichletBasal(Enth0); CHKERRQ(ierr);
        } else if (base_is_cold) {
          // cold, grounded base case:  Neumann q . n = q_lith . n + F_b   and   q = - K_i \nabla H
          ierr = esys.setNeumannBasal(- (G0(i,j) + (*Rb)(i,j)) / ice_K); CHKERRQ(ierr);
        } else {
          // warm, grounded base case
          if (k1_istemperate) {
            // positive thickness of temperate ice:  Neumann q . n = 0 and q = - K_0 \nabla H
            //   so H(k=1)-H(k=0) = 0
            ierr = esys.setNeumannBasal(0.0); CHKERRQ(ierr);
          } else {
            // no thickness of temperate ice:  Dirichlet  H = H_s(pbasal)
            ierr = esys.setDirichletBasal(esys.Enth_s[0]); CHKERRQ(ierr);
          }
        }

        // solve the system
        PetscErrorCode pivoterr;
        ierr = esys.solveThisColumn(&Enthnew,pivoterr); CHKERRQ(ierr);
        if (pivoterr != 0) {
          ierr = PetscPrintf(PETSC_COMM_SELF,
            "\n\ntridiagonal solve of enthSystemCtx in enthalpyAndDrainageStep() FAILED at (%d,%d)\n"
                " with zero pivot position %d; viewing system to m-file ... \n",
            i, j, pivoterr); CHKERRQ(ierr);
          ierr = esys.reportColumnZeroPivotErrorMFile(pivoterr); CHKERRQ(ierr);
          SETERRQ(1,"PISM ERROR in enthalpyDrainageStep()\n");
        }
        if (viewOneColumn && issounding(i,j)) {
          ierr = PetscPrintf(PETSC_COMM_SELF,
            "\n\nin enthalpyAndDrainageStep(): viewing enthSystemCtx at (i,j)=(%d,%d) to m-file ... \n\n",
            i, j); CHKERRQ(ierr);
          ierr = esys.viewColumnInfoMFile(Enthnew, fMz); CHKERRQ(ierr);
        }

        // thermodynamic basal melt rate causes water to be added to layer
        PetscScalar bwatnew = vbwat(i,j);
        if (mask.grounded(i,j)) {
          bwatnew += vbmr(i,j) * dt_secs;
        }

        // drain ice segments by mechanism in [\ref AschwandenBuelerKhroulevBlatter],
        //   using DrainageCalculator dc
        PetscScalar Hdrainedtotal = 0.0;
        for (PetscInt k=0; k < ks; k++) {
          if (Enthnew[k] > esys.Enth_s[k]) { // avoid doing any more work if cold
            if (Enthnew[k] >= esys.Enth_s[k] + 0.5 * L) {
              liquifiedCount++; // count these rare events ...
              Enthnew[k] = esys.Enth_s[k] + 0.5 * L; //  but lose the energy
            }
            const PetscReal p = EC->getPressureFromDepth(vH(i,j) - fzlev[k]); // FIXME task #7297
            PetscReal omega;
            EC->getWaterFraction(Enthnew[k], p, omega);  // return code not checked
            if (omega > 0.01) {
              PetscReal fractiondrained = dc.get_drainage_rate(omega) * dt_secs; // pure number
              fractiondrained = PetscMin(fractiondrained, omega - 0.01); // only drain down to 0.01
              Hdrainedtotal += fractiondrained * fdz;  // always a positive contribution
              Enthnew[k] -= fractiondrained * L;
            }
          }
        }

        // in grounded case, add to both basal melt rate and bwat; if floating,
        // Hdrainedtotal is discarded because ocean determines basal melt rate
        if (mask.grounded(i,j)) {
          vbmr(i,j) += Hdrainedtotal / dt_secs;
          bwatnew += Hdrainedtotal;
        }

        // finalize Enthnew[]:  apply bulge limiter and transfer column
        //   into vWork3d; communication will occur later
        const PetscReal lowerEnthLimit = Enth_ks - bulgeEnthMax;
        for (PetscInt k=0; k < ks; k++) {
          if (Enthnew[k] < lowerEnthLimit) {
            *bulgeCount += 1;      // count the columns which have very large cold 
            Enthnew[k] = lowerEnthLimit;  // limit advection bulge ... enthalpy not too low
          }
        }
        ierr = vWork3d.setValColumnPL(i,j,Enthnew); CHKERRQ(ierr);

        // finalize bwat value
        bwatnew -= bwat_decay_rate * dt_secs;
        if (is_floating) {
          // if floating assume maximally saturated till to avoid "shock" if grounding line advances
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS/GAIN (TO/FROM OCEAN)!!
          vbwat(i,j) = bwat_max;
        } else {
          // limit bwat to be in [0.0, bwat_max]
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO INFINITY AND BEYOND)!!
          vbwat(i,j) = PetscMax(0.0, PetscMin(bwat_max, bwatnew) );
        }

      } // end explicit scoping
      
      donewithcolumn: 
      { }  // odd thing: something needs to follow goto target to get compilation

    }
  }

  ierr = artm.end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.end_access(); CHKERRQ(ierr);

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vbwat.end_access(); CHKERRQ(ierr);
  ierr = Rb->end_access(); CHKERRQ(ierr);
  ierr = G0.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);
  ierr = liqfrac_surface.end_access(); CHKERRQ(ierr);

  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = Sigma3->end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vWork3d.end_access(); CHKERRQ(ierr);

  delete [] Enthnew;

  *liquifiedVol = ((double) liquifiedCount) * fdz * grid.dx * grid.dy;
  return 0;
}



#include "coupler/PCFactory.hh"
#include "coupler/atmosphere/PISMAtmosphere.hh"
#include "coupler/surface/PISMSurface.hh"
#include "coupler/ocean/PISMOcean.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "VARKPISMR %s (modified basic evolution run mode: variable cold-ice conductivity)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    ierr = check_old_option_and_stop(com, "-boot_from", "-boot_file"); CHKERRQ(ierr); 

    bool iset, bfset;
    ierr = PISMOptionsIsSet("-i", iset); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", bfset); CHKERRQ(ierr);
    string usage =
      "  varkpismr {-i IN.nc|-boot_file IN.nc} [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "notes:\n"
      "  * one of -i or -boot_file is required\n"
      "  * if -boot_file is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
      ierr = PetscPrintf(com,
         "\nPISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "varkpismr", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "varkpismr", required, usage.c_str()); CHKERRQ(ierr);
    }

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid g(com, rank, size, config);
    IcevarkModel m(g, config, overrides);

    // Initialize boundary models:
    PAFactory pa(g, config);
    PISMAtmosphereModel *atmosphere;

    PSFactory ps(g, config);
    PISMSurfaceModel *surface;

    POFactory po(g, config);
    PISMOceanModel *ocean;

    ierr = PetscOptionsBegin(com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);
    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    surface->attach_atmosphere_model(atmosphere);

    m.attach_ocean_model(ocean);
    m.attach_surface_model(surface);
    ierr = m.setExecName("varkpismr"); CHKERRQ(ierr);

    ierr = m.init(); CHKERRQ(ierr);

    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);
    // provide a default output file name if no -o option is given.
    ierr = m.writeFiles("varkunnamed.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
