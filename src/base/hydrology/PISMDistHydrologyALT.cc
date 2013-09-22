// Copyright (C) 2012-2013 PISM Authors
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

#include "PISMHydrology.hh"
#include "PISMVars.hh"
#include "pism_options.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"


//! Update the model state variables W,P by running the subglacial hydrology model, but use a different update method for pressure than is in PISMDistributedHydrology::update().
/*! In contrast to PISMDistributedHydrology::update() this method uses the
hydraulic potential (i.e. it computes the gradient of it) but it does not use
Qstag, the staggered-grid values of the advective flux. */
PetscErrorCode PISMDistHydrologyALT::update(PetscReal icet, PetscReal icedt) {
  PetscErrorCode ierr;

  // if asked for the identical time interval versus last time, then
  //   do nothing; otherwise assume that [my_t,my_t+my_dt] is the time
  //   interval on which we are solving
  if ((fabs(icet - t) < 1e-12) && (fabs(icedt - dt) < 1e-12))
    return 0;
  // update PISMComponent times: t = current time, t+dt = target time
  t = icet;
  dt = icedt;

  // make sure W,P have valid ghosts before starting hydrology steps
  ierr = W.update_ghosts(); CHKERRQ(ierr);
  ierr = P.update_ghosts(); CHKERRQ(ierr);

  // from current ice geometry/velocity variables, initialize Po and cbase
  ierr = update_cbase(cbase); CHKERRQ(ierr);

  const PetscReal
            rg = config.get("fresh_water_density") * config.get("standard_gravity"),
            nglen = config.get("Glen_exponent"),
            Aglen = config.get("ice_softness"),
            c1 = config.get("hydrology_cavitation_opening_coefficient"),
            c2 = config.get("hydrology_creep_closure_coefficient"),
            Wr = config.get("hydrology_roughness_scale"),
            phi0 = config.get("hydrology_regularizing_porosity");

  const PetscReal  omegax = 1.0 / (grid.dx * grid.dx),
                   omegay = 1.0 / (grid.dy * grid.dy);

  PetscReal ht = t, hdt, // hydrology model time and time step
            maxKW, maxV, maxD;
  PetscReal icefreelost = 0, oceanlost = 0, negativegain = 0, nullstriplost = 0,
            delta_icefree, delta_ocean, delta_neggain, delta_nullstrip;

  PetscReal PtoCFLratio,  // for reporting ratio of dtCFL to dtDIFFP
            cumratio = 0.0;
  PetscInt hydrocount = 0; // count hydrology time steps

  while (ht < t + dt) {
    hydrocount++;

#if (PISM_DEBUG==1)
    ierr = check_water_thickness_nonnegative(W); CHKERRQ(ierr);
    ierr = check_Wtil_bounds(); CHKERRQ(ierr);
#endif

    // update Wtilnew (the actual step) from W and Wtil
    ierr = raw_update_Wtil(hdt); CHKERRQ(ierr);
    ierr = boundary_mass_changes(Wtilnew, delta_icefree, delta_ocean,
                                 delta_neggain, delta_nullstrip); CHKERRQ(ierr);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // note that ice dynamics can change overburden pressure, so we can only check P
    //   bounds if thk has not changed; if thk could have just changed, such as in the
    //   first time through the current loop, we enforce them
    ierr = check_P_bounds((hydrocount == 1)); CHKERRQ(ierr);

    ierr = subglacial_hydraulic_potential(psi); CHKERRQ(ierr);
    ierr = psi.update_ghosts(); CHKERRQ(ierr);

    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);
    ierr = Wstag.update_ghosts(); CHKERRQ(ierr);

    ierr = conductivity_staggered(Kstag,maxKW); CHKERRQ(ierr);
    ierr = Kstag.update_ghosts(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);

    ierr = adaptive_for_WandP_evolution(ht, t+dt, maxKW, hdt, maxV, maxD, PtoCFLratio); CHKERRQ(ierr);
    cumratio += PtoCFLratio;

    if ((inputtobed != NULL) || (hydrocount==1)) {
      ierr = get_input_rate(ht,hdt,total_input); CHKERRQ(ierr);
    }

    // update Pnew from time step
    const PetscReal  CC = (rg * hdt) / phi0;
    PetscReal  Open, Close, divflux, ZZ,
               dpsie, dpsiw, dpsin, dpsis;
    ierr = overburden_pressure(Pover); CHKERRQ(ierr);

    MaskQuery M(*mask);

    ierr = P.begin_access(); CHKERRQ(ierr);
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wtil.begin_access(); CHKERRQ(ierr);
    ierr = Wtilnew.begin_access(); CHKERRQ(ierr);
    ierr = cbase.begin_access(); CHKERRQ(ierr);
    ierr = psi.begin_access(); CHKERRQ(ierr);
    ierr = Wstag.begin_access(); CHKERRQ(ierr);
    ierr = Kstag.begin_access(); CHKERRQ(ierr);
    ierr = total_input.begin_access(); CHKERRQ(ierr);
    ierr = mask->begin_access(); CHKERRQ(ierr);
    ierr = Pover.begin_access(); CHKERRQ(ierr);
    ierr = Pnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (M.ice_free_land(i,j))
          Pnew(i,j) = 0.0;
        else if (M.ocean(i,j))
          Pnew(i,j) = Pover(i,j);
        else if (W(i,j) <= 0.0)
          Pnew(i,j) = Pover(i,j);
        else {
          // opening and closure terms in pressure equation
          Open = PetscMax(0.0,c1 * cbase(i,j) * (Wr - W(i,j)));
          Close = c2 * Aglen * pow(Pover(i,j) - P(i,j),nglen) * W(i,j);

// ALT start
          // divergence of flux
          const bool knowne = (M.ice_free_land(i+1,j) || M.ocean(i+1,j)),
                     knownw = (M.ice_free_land(i-1,j) || M.ocean(i-1,j)),
                     knownn = (M.ice_free_land(i,j+1) || M.ocean(i,j+1)),
                     knowns = (M.ice_free_land(i,j-1) || M.ocean(i,j-1));
          dpsie = psi(i+1,j) - psi(i,j);
          dpsiw = psi(i,j)   - psi(i-1,j);
          dpsin = psi(i,j+1) - psi(i,j);
          dpsis = psi(i,j)   - psi(i,j-1);
          if (stripwidth > 0.0) {
            const bool nullij = (in_null_strip(i,j));
            if (nullij || in_null_strip(i+1,j))
              dpsie = 0.0;
            if (nullij || in_null_strip(i-1,j))
              dpsiw = 0.0;
            if (nullij || in_null_strip(i,j+1))
              dpsin = 0.0;
            if (nullij || in_null_strip(i,j-1))
              dpsis = 0.0;
          }
          divflux = 0.0;
          if (!knowne && !knownw) {
            const PetscReal We = Wstag(i,  j,0),
                            Ww = Wstag(i-1,j,0),
                            Ke = Kstag(i,  j,0),
                            Kw = Kstag(i-1,j,0);
            divflux += omegax * ( Ke * We * dpsie - Kw * Ww * dpsiw );
          }
          if (!knownn && !knowns) {
            const PetscReal Wn = Wstag(i,j  ,1),
                            Ws = Wstag(i,j-1,1);
            const PetscReal Kn = Kstag(i,j  ,1),
                            Ks = Kstag(i,j-1,1);
            divflux += omegay * ( Kn * Wn * dpsin - Ks * Ws * dpsis );
          }
// ALT end
          // pressure update equation
          ZZ = Close - Open + total_input(i,j) - (Wtilnew(i,j) - Wtil(i,j)) / hdt;
          Pnew(i,j) = P(i,j) + CC * ( divflux + ZZ );
          // projection to enforce  0 <= P <= P_o
          Pnew(i,j) = PetscMin(PetscMax(0.0, Pnew(i,j)), Pover(i,j));
        }
      }
    }
    ierr = P.end_access(); CHKERRQ(ierr);
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wtil.end_access(); CHKERRQ(ierr);
    ierr = Wtilnew.end_access(); CHKERRQ(ierr);
    ierr = cbase.end_access(); CHKERRQ(ierr);
    ierr = Pnew.end_access(); CHKERRQ(ierr);
    ierr = Pover.end_access(); CHKERRQ(ierr);
    ierr = psi.end_access(); CHKERRQ(ierr);
    ierr = total_input.end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = Kstag.end_access(); CHKERRQ(ierr);
    ierr = mask->end_access(); CHKERRQ(ierr);

    // FIXME: following chunk is code duplication with PISMRoutingHydrology::update()

    // update Wnew (the actual step) from W, Wtil, Wtilnew, Wstag, Qstag, total_input
    ierr = raw_update_W(hdt); CHKERRQ(ierr);
    ierr = boundary_mass_changes(Wnew, delta_icefree, delta_ocean,
                                 delta_neggain, delta_nullstrip); CHKERRQ(ierr);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // transfer new into old
    ierr = Wnew.update_ghosts(W); CHKERRQ(ierr);
    ierr = Wtilnew.copy_to(Wtil); CHKERRQ(ierr);
    ierr = Pnew.update_ghosts(P); CHKERRQ(ierr);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  if (report_mass_accounting) {
    ierr = verbPrintf(2, grid.com,
                      " 'distributed' hydrology summary:\n"
                      "     %d hydrology sub-steps with average dt = %.7f years = %.2f s\n"
                      "        (average of %.2f steps per CFL time; last max |V| = %.2e m s-1; last max D = %.2e m^2 s-1)\n"
                      "     ice free land lost = %.3e kg, ocean lost = %.3e kg\n"
                      "     negative bmelt gain = %.3e kg, null strip lost = %.3e kg\n",
                      hydrocount, grid.convert(dt/hydrocount, "seconds", "years"), dt/hydrocount,
                      cumratio/hydrocount, maxV, maxD,
                      icefreelost, oceanlost,
                      negativegain, nullstriplost); CHKERRQ(ierr);
  }
  return 0;
}

