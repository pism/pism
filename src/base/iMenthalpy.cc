// Copyright (C) 2009-2013 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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

#include "iceModel.hh"
#include "enthSystem.hh"
#include "DrainageCalculator.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"
#include "PISMHydrology.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "bedrockThermalUnit.hh"
#include "enthalpyConverter.hh"
#include "pism_options.hh"

//! \file iMenthalpy.cc Methods of IceModel which implement the enthalpy formulation of conservation of energy.


//! Compute Enth3 from temperature T3 by assuming the ice has zero liquid fraction.
/*!
First this method makes sure the temperatures is at most the pressure-melting
value, before computing the enthalpy for that temperature, using zero liquid
fraction.

Because of how EnthalpyConverter::getPressureFromDepth() works, the energy
content in the air is set to the value that ice would have if it a chunk of it
occupied the air; the atmosphere actually has much lower energy content.  It is
done this way for regularity (i.e. dEnth/dz computations).

Because Enth3 gets set, does ghost communication to finish.
 */
PetscErrorCode IceModel::compute_enthalpy_cold(IceModelVec3 &temperature, IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = temperature.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  PetscScalar *Tij, *Enthij; // columns of these values
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = temperature.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = result.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME issue #15
        ierr = EC->getEnthPermissive(Tij[k],0.0,EC->getPressureFromDepth(depth),
                                    Enthij[k]); CHKERRQ(ierr);
      }
    }
  }

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = temperature.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.update_ghosts(); CHKERRQ(ierr);

  return 0;
}


//! Compute Enth3 from temperature T3 and liquid fraction.
/*!
Because Enth3 gets set, does ghost communication to finish.
 */
PetscErrorCode IceModel::compute_enthalpy(IceModelVec3 &temperature,
                                          IceModelVec3 &liquid_water_fraction,
                                          IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = temperature.begin_access(); CHKERRQ(ierr);
  ierr = liquid_water_fraction.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  PetscScalar *Tij, *Liqfracij, *Enthij; // columns of these values
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = temperature.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = liquid_water_fraction.getInternalColumn(i,j,&Liqfracij); CHKERRQ(ierr);
      ierr = result.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME issue #15
        ierr = EC->getEnthPermissive(Tij[k],Liqfracij[k],
                      EC->getPressureFromDepth(depth), Enthij[k]); CHKERRQ(ierr);
      }
    }
  }

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = temperature.end_access(); CHKERRQ(ierr);
  ierr = liquid_water_fraction.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

//! Compute the liquid fraction corresponding to Enth3, and put in a global IceModelVec3 provided by user.
/*!
Does not communicate ghosts for IceModelVec3 result
 */
PetscErrorCode IceModel::compute_liquid_water_fraction(IceModelVec3 &enthalpy,
                                                       IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = result.set_name("liqfrac"); CHKERRQ(ierr);
  ierr = result.set_attrs(
     "diagnostic",
     "liquid water fraction in ice (between 0 and 1)",
     "1", ""); CHKERRQ(ierr);

  PetscScalar *omegaij, *Enthij; // columns of these values
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = enthalpy.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result.getInternalColumn(i,j,&omegaij); CHKERRQ(ierr);
      ierr = enthalpy.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME issue #15
        ierr = EC->getWaterFraction(Enthij[k],EC->getPressureFromDepth(depth),
                                   omegaij[k]); CHKERRQ(ierr);
      }
    }
  }
  ierr = enthalpy.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  return 0;
}

//! Compute the CTS field, CTS = E/E_s(p), from Enth3, and put in a global IceModelVec3 provided by user.
/*!
The actual cold-temperate transition surface (CTS) is the level set CTS = 1.

Does not communicate ghosts for IceModelVec3 result.
 */
PetscErrorCode IceModel::setCTSFromEnthalpy(IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = result.set_name("cts"); CHKERRQ(ierr);
  ierr = result.set_attrs(
     "diagnostic",
     "cts = E/E_s(p), so cold-temperate transition surface is at cts = 1",
     "", ""); CHKERRQ(ierr);

  PetscScalar *CTSij, *Enthij; // columns of these values
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result.getInternalColumn(i,j,&CTSij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME issue #15
        CTSij[k] = EC->getCTS(Enthij[k], EC->getPressureFromDepth(depth));
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  return 0;
}

//! Update ice enthalpy field based on conservation of energy.
/*!
This method is documented by the page \ref bombproofenth and by [\ref
AschwandenBuelerKhroulevBlatter].

This method updates IceModelVec3 vWork3d = vEnthnew and IceModelVec2S vbmr.
No communication of ghosts is done for any of these fields.

We use an instance of enthSystemCtx.

Regarding drainage, see [\ref AschwandenBuelerKhroulevBlatter] and references therein.

\image html BC-decision-chart.png "Setting the basal boundary condition"
 */
PetscErrorCode IceModel::enthalpyAndDrainageStep(PetscScalar* vertSacrCount,
                                                 PetscScalar* liquifiedVol,
                                                 PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  assert(config.get_flag("do_cold_ice_methods") == false);

  // essentially physical constants:
  const PetscScalar
    ice_rho      = config.get("ice_density"), // kg m-3
    L            = config.get("water_latent_heat_fusion"), // J kg-1
    // constants controlling the numerical method:
    bulgeEnthMax = config.get("enthalpy_cold_bulge_max"); // J kg-1

  bool viewOneColumn;
  ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

  DrainageCalculator dc(config);

  IceModelVec2S *Rb;
  IceModelVec3 *u3, *v3, *w3, *strain_heating3;
  ierr = stress_balance->get_basal_frictional_heating(Rb); CHKERRQ(ierr);
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);
  ierr = stress_balance->get_volumetric_strain_heating(strain_heating3); CHKERRQ(ierr);

  PetscScalar *Enthnew;
  Enthnew = new PetscScalar[grid.Mz_fine];  // new enthalpy in column

  enthSystemCtx esys(config, Enth3, grid.dx, grid.dy, dt_TempAge,
                     grid.dz_fine, grid.Mz_fine, "enth", EC);

  // Now get map-plane coupler fields: Dirichlet upper surface
  // boundary and mass balance lower boundary under shelves
  assert(surface != NULL);
  ierr = surface->ice_surface_temperature(ice_surface_temp); CHKERRQ(ierr);
  ierr = surface->ice_surface_liquid_water_fraction(liqfrac_surface); CHKERRQ(ierr);

  assert(ocean != NULL);
  ierr = ocean->shelf_base_mass_flux(shelfbmassflux); CHKERRQ(ierr);
  ierr = ocean->shelf_base_temperature(shelfbtemp); CHKERRQ(ierr);

  IceModelVec2S basal_heat_flux = vWork2d[0];
  ierr = basal_heat_flux.set_attrs("internal", "upward heat flux at z=0",
                                   "W m-2", ""); CHKERRQ(ierr);
  assert(btu != NULL);
  ierr = btu->get_upward_geothermal_flux(basal_heat_flux); CHKERRQ(ierr);

  IceModelVec2S basal_water_thickness = vWork2d[1];
  ierr = basal_water_thickness.set_attrs("internal", "current amount of basal water",
                                         "m", ""); CHKERRQ(ierr);
  assert(subglacial_hydrology != NULL);
  ierr = subglacial_hydrology->subglacial_water_thickness(basal_water_thickness); CHKERRQ(ierr);

  ierr = ice_surface_temp.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

  // get other map-plane fields
  ierr = liqfrac_surface.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  ierr = Rb->begin_access(); CHKERRQ(ierr);
  ierr = basal_heat_flux.begin_access(); CHKERRQ(ierr);
  ierr = basal_water_thickness.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  // these are accessed a column at a time
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = strain_heating3->begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vWork3d.begin_access(); CHKERRQ(ierr);

  PetscInt liquifiedCount = 0;

  MaskQuery mask(vMask);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // ignore advection and strain heating in ice if isMarginal
      const bool isMarginal = checkThinNeigh(vH(i+1, j), vH(i+1, j+1), vH(i, j+1), vH(i-1, j+1),
                                             vH(i-1, j), vH(i-1, j-1), vH(i, j-1), vH(i+1, j-1));

      ierr = esys.initThisColumn(i, j, isMarginal, vH(i, j),
                                  u3, v3, w3, strain_heating3); CHKERRQ(ierr);

      const int ks = esys.ks();
      const bool
        ice_free_column = (ks == 0),
        is_floating     = mask.ocean(i, j);

      // enthalpy and pressures at top of ice
      const PetscScalar
        depth_ks = vH(i, j) - ks * grid.dz_fine,
        p_ks     = EC->getPressureFromDepth(depth_ks); // FIXME issue #15

      PetscScalar Enth_ks;
      ierr = EC->getEnthPermissive(ice_surface_temp(i, j), liqfrac_surface(i, j),
                                   p_ks, Enth_ks); CHKERRQ(ierr);

      // deal completely with columns with no ice; enthalpy and vbmr need setting
      if (ice_free_column) {
        ierr = vWork3d.setColumn(i, j, Enth_ks); CHKERRQ(ierr);
        if (mask.floating_ice(i, j)) {
          vbmr(i, j) = shelfbmassflux(i, j);
        } else {
          vbmr(i, j) = 0.0;    // no basal melt rate on ice free land
        }
        continue;
      } // end of if (ice_free_column)

      if (esys.lambda() < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1

      bool base_is_cold = (esys.Enth[0] < esys.Enth_s[0]);
      const PetscScalar p1 = EC->getPressureFromDepth(vH(i, j) - grid.dz_fine); // FIXME issue #15
      const bool k1_istemperate = EC->isTemperate(esys.Enth[1], p1); // level  z = + \Delta z

      // if there is subglacial water, don't allow ice base enthalpy to be below
      // pressure-melting; that is, assume subglacial water is at the pressure-
      // melting temperature and enforce continuity of temperature
      if (basal_water_thickness(i, j) > 0.0 && base_is_cold) {
        esys.Enth[0] = esys.Enth_s[0];
        base_is_cold = false;
      }

      // can now determine melt, but only preliminarily because of drainage,
      //   from heat flux out of bedrock, heat flux into ice, and frictional heating
      if (is_floating) {
        vbmr(i, j) = shelfbmassflux(i, j);
      } else {
        if (base_is_cold) {
          vbmr(i, j) = 0.0;  // zero melt rate if cold base
        } else {
          const PetscScalar p_0 = EC->getPressureFromDepth(vH(i, j)); // FIXME issue #15
          PetscScalar hf_up;
          if (k1_istemperate) {
            const PetscScalar
              Tpmp_0 = EC->getMeltingTemp(p_0),
              Tpmp_1 = EC->getMeltingTemp(p1);

            hf_up = -esys.k_from_T(Tpmp_0) * (Tpmp_1 - Tpmp_0) / grid.dz_fine;
          } else {
            PetscScalar T_0;
            ierr = EC->getAbsTemp(esys.Enth[0], p_0, T_0); CHKERRQ(ierr);
            const PetscScalar K_0 = esys.k_from_T(T_0) / EC->c_from_T(T_0);

            hf_up = - K_0 * (esys.Enth[1] - esys.Enth[0]) / grid.dz_fine;
          }

          // compute basal melt rate from flux balance:
          //
          // vbmr = - Mb / rho in [\ref AschwandenBuelerKhroulevBlatter];
          //
          // after we compute it we make sure there is no refreeze if
          // there is no available basal water
          vbmr(i, j) = ( (*Rb)(i, j) + basal_heat_flux(i, j) - hf_up ) / (ice_rho * L);

          if (basal_water_thickness(i, j) <= 0 && vbmr(i, j) < 0)
            vbmr(i, j) = 0.0;
        }
      }

      ierr = esys.setDirichletSurface(Enth_ks); CHKERRQ(ierr);

      // determine lowest-level equation at bottom of ice; see decision chart
      //   in [\ref AschwandenBuelerKhroulevBlatter], and page documenting BOMBPROOF
      if (is_floating) {
        // floating base: Dirichlet application of known temperature from ocean
        //   coupler; assumes base of ice shelf has zero liquid fraction
        PetscScalar Enth0;
        ierr = EC->getEnthPermissive(shelfbtemp(i, j), 0.0, EC->getPressureFromDepth(vH(i, j)),
                                     Enth0); CHKERRQ(ierr);
        ierr = esys.setDirichletBasal(Enth0); CHKERRQ(ierr);
      } else if (base_is_cold) {
        // cold, grounded base (Neumann) case:  q . n = q_lith . n + F_b
        ierr = esys.setBasalHeatFlux(basal_heat_flux(i, j) + (*Rb)(i, j)); CHKERRQ(ierr);
      } else {
        // warm, grounded base case
        ierr = esys.setBasalHeatFlux(0.0); CHKERRQ(ierr);
      }

      // solve the system
      ierr = esys.solveThisColumn(&Enthnew); CHKERRQ(ierr);

      if (viewOneColumn && issounding(i, j)) {
        ierr = esys.viewColumnInfoMFile(Enthnew, grid.Mz_fine); CHKERRQ(ierr);
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
          const PetscReal depth = vH(i, j) - k * grid.dz_fine,
            p = EC->getPressureFromDepth(depth); // FIXME issue #15
          PetscReal omega;
          EC->getWaterFraction(Enthnew[k], p, omega);  // return code not checked
          if (omega > 0.01) {
            PetscReal fractiondrained = dc.get_drainage_rate(omega) * dt_TempAge; // pure number
            fractiondrained = PetscMin(fractiondrained, omega - 0.01); // only drain down to 0.01
            Hdrainedtotal += fractiondrained * grid.dz_fine;  // always a positive contribution
            Enthnew[k] -= fractiondrained * L;
          }
        }
      }

      // in grounded case, add to basal melt rate; if floating,
      // Hdrainedtotal is discarded because ocean determines basal melt rate
      if (mask.grounded(i, j)) {
        vbmr(i, j) += Hdrainedtotal / dt_TempAge;
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
      ierr = vWork3d.setValColumnPL(i, j, Enthnew); CHKERRQ(ierr);

    }
  }

  ierr = ice_surface_temp.end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.end_access(); CHKERRQ(ierr);

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = Rb->end_access(); CHKERRQ(ierr);
  ierr = basal_heat_flux.end_access(); CHKERRQ(ierr);
  ierr = basal_water_thickness.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);
  ierr = liqfrac_surface.end_access(); CHKERRQ(ierr);

  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = strain_heating3->end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vWork3d.end_access(); CHKERRQ(ierr);

  delete [] Enthnew;

  *liquifiedVol = ((double) liquifiedCount) * grid.dz_fine * grid.dx * grid.dy;
  return 0;
}

/*
  The decision above was produced using this TikZ source code:

% Define block styles
\tikzstyle{decision} = [ellipse, draw, text width=7em, text badly centered, inner sep=2pt]
\tikzstyle{block} = [rectangle, draw, text width=5em, text badly centered, rounded corners, minimum height=4em]
\tikzstyle{line} = [draw, -latex']

\begin{tikzpicture}[node distance = 3cm, auto]
    % Place nodes
    \node (invisiblestart) {};

    \node [decision, below of=invisiblestart, text height=0.2cm] (coldvstemp) {$H<H_{\text s}(p)$ ?};
    \node [decision, left of=coldvstemp, xshift=-4em] (excludebad) {$\eta_{\text b}>0$ ?};
    \node [block, below of=excludebad, text width=6em] (fixbad) {$H := H_{\text s}(p)$};

    % edges
    \path [line] (invisiblestart) -- (coldvstemp);
    \path [line] (excludebad) -- node [text width=6em] {yes (consider base to be temperate)} (fixbad);

    % cold branch:
    \node [block, left of=fixbad, text width=7.5em] (coldmodeltype) {Eqn (49) is Neumann b.c.~for Eqn (67); $M_b=0$};
    % edges
    \path [line] (coldvstemp) -- node {yes} (excludebad);
    \path [line] (excludebad) -- node {no} (coldmodeltype);

    % temperate branch
    \node [block, below of=coldvstemp, text width=12em] (qtemperate) {$\nabla H \cdot \bn=0$ is Neumann b.c.~for Eqn (67)};
    \node [decision, below left of=qtemperate, text width=8em] (tempthick) {positive thickness of temperate ice at base?};
    \node [block, below right of=tempthick, text width=10em] (Mbforqtemperate) {$\bq = - k(H,p)\nabla T_{\text m}(p)$ \\ at ice base};
    \node [block, below left of=tempthick, text width=9em, xshift=-4em] (Mbforqcold) {$\bq = - K_{\text i}(H) \nabla H$ \\ at ice base};
    \node [block, below left of=Mbforqtemperate, text width=9em] (getMbtemp) {compute $M_b$ from Eqn (50) or (66)};

    % edges
    \path [line] (fixbad) -- (qtemperate);
    \path [line] (coldvstemp) -- node {no} (qtemperate);
    \path [line] (tempthick) -- node [above] {no} (Mbforqcold);
    \path [line] (tempthick) -- node {yes} (Mbforqtemperate);
    \path [line] (qtemperate) -- (tempthick);
    \path [line] (Mbforqcold) -- node {} (getMbtemp);
    \path [line] (Mbforqtemperate) -- node {} (getMbtemp);
\end{tikzpicture}
 */
