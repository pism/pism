// Copyright (C) 2004--2013 Torsten Albrecht and Constantine Khroulev
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
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

#include <cmath>
#include <petscdmda.h>

#include "iceModel.hh"
#include "Mask.hh"
#include "PISMOcean.hh"

//! \file iMicebergs.cc Methods implementing PIK option -kill_icebergs [\ref Winkelmannetal2011].


//! \brief Identify and eliminate free-floating icebergs, which cause
//! well-posedness (invertibility) problems for stress solvers.
/*!

 * Icebergs are, in this context, floating regions that are \e not attached,
 * through a chain of positive thickness ice-filled cells, to at least one
 * grounded cell. They are observed to cause unrealistically large velocities
 * that may (numerically) affect the ice velocities everywhere. They cause the
 * SSA operator to have a nontrivial null space, or, under approximation
 * errors, they lead to extremely small time steps and can eventually cause a
 * KSP-ERROR.

 * This method calls the routines which identify and then eliminate these icebergs.

 * FIXME: a fundamental aspect of the semantics here is not clear to me
 * (bueler), namely how many times the iceberg-eliminate "sweep" might occur,
 * and what properties control that? for now, should there be some
 * (low-verbosity) indication that it is occurring, such as when more than one
 * sweep happened?

 * FIXME: this package of methods *might* appropriately be a class

 * FIXME: this package of routines *should* have a regression test
 */
PetscErrorCode IceModel::killIceBergs() {
  PetscErrorCode ierr;

  ierr = findIceBergCandidates(); CHKERRQ(ierr);
  ierr = identifyNotAnIceBerg(); CHKERRQ(ierr);
  ierr = killIdentifiedIceBergs(); CHKERRQ(ierr);
  ierr = killEasyIceBergs(); CHKERRQ(ierr);
  return 0;
}


/*!
 * The aim of this routine is to find floating regions that *might* be
 * icebergs. If these regions actually *are* icebergs is checked in
 * identifyNotAnIceBerg().
 */
PetscErrorCode IceModel::findIceBergCandidates() {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "######### findIceBergCandidates() start\n"); CHKERRQ(ierr);

  const PetscInt Mx = grid.Mx, My = grid.My;

  PetscReal sea_level;
  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else { SETERRQ(grid.com, 2, "PISM ERROR: ocean == NULL"); }

  double ocean_rho = config.get("sea_water_density"),
    ice_rho = config.get("ice_density");

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      const PetscScalar hgrounded = vbed(i, j) + vH(i, j),
        hfloating = sea_level + (1.0 - ice_rho / ocean_rho) * vH(i, j);

      //cut of border of computational domain
      if (hgrounded < hfloating && (i <= 0 || i >= Mx - 1 || j <= 0 || j >= My - 1)) {
        vH(i, j) = 0.0;
        vIcebergMask(i, j) = ICEBERGMASK_STOP_OCEAN;
        vMask(i, j) = MASK_ICE_FREE_OCEAN;
      } else {
        vIcebergMask(i, j) = ICEBERGMASK_NOT_SET;
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);

  ierr = vMask.update_ghosts(); CHKERRQ(ierr);
  ierr = vIcebergMask.update_ghosts(); CHKERRQ(ierr);

  // set all floating points to ICEBERGMASK_ICEBERG_CAND
  MaskQuery M(vMask);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      if (vIcebergMask(i, j) == ICEBERGMASK_NOT_SET) {
        if (M.floating_ice(i, j)) {
          vIcebergMask(i, j) = ICEBERGMASK_ICEBERG_CAND;
        }
      }
    }
  }
  ierr = vIcebergMask.end_access(); CHKERRQ(ierr);

  ierr = vIcebergMask.update_ghosts(); CHKERRQ(ierr);

  // set borders of shelves/icebergs to ICEBERGMASK_STOP_ATTACHED or ICEBERGMASK_STOP_OCEAN respectively.
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      if (vIcebergMask(i, j) == ICEBERGMASK_NOT_SET) {

        planeStar<int> mask = vIcebergMask.int_star(i, j);

        bool neighbor_is_candidate = ( mask.e == ICEBERGMASK_ICEBERG_CAND ||
                                       mask.w == ICEBERGMASK_ICEBERG_CAND ||
                                       mask.n == ICEBERGMASK_ICEBERG_CAND ||
                                       mask.s == ICEBERGMASK_ICEBERG_CAND ||
                                       // checks below should not be needed
                                       vIcebergMask(i + 1, j + 1) == ICEBERGMASK_ICEBERG_CAND ||
                                       vIcebergMask(i + 1, j - 1) == ICEBERGMASK_ICEBERG_CAND ||
                                       vIcebergMask(i - 1, j + 1) == ICEBERGMASK_ICEBERG_CAND ||
                                       vIcebergMask(i - 1, j - 1) == ICEBERGMASK_ICEBERG_CAND);

        if (M.grounded(i, j) && neighbor_is_candidate)
          vIcebergMask(i, j) = ICEBERGMASK_STOP_ATTACHED;
        else if (M.ice_free_ocean(i, j) && neighbor_is_candidate)
          vIcebergMask(i, j) = ICEBERGMASK_STOP_OCEAN;
      }
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.end_access(); CHKERRQ(ierr);

  ierr = vIcebergMask.update_ghosts(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::identifyNotAnIceBerg() {
  PetscErrorCode ierr;

  ierr = verbPrintf(4, grid.com, "######### identifyNotAnIceBerg() start\n"); CHKERRQ(ierr);

  // this communication of ghost values is done here to make sure that asking
  // about neighbouring values in this routine doesn't lead to inconsistencies
  // in parallel computation, if the neighbour belongs to another processor
  // domain.
  // FIXME: this is probably redundant
  ierr = vIcebergMask.update_ghosts(); CHKERRQ(ierr);

  bool done = false;
  PetscInt loopcount = 0;
  while(! done){

    done = true;

    ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);
    for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

        planeStar<int> mask = vIcebergMask.int_star(i, j);

        bool attached_to_grounded = (mask.e == ICEBERGMASK_STOP_ATTACHED ||
                                     mask.w == ICEBERGMASK_STOP_ATTACHED ||
                                     mask.n == ICEBERGMASK_STOP_ATTACHED ||
                                     mask.s == ICEBERGMASK_STOP_ATTACHED),

          attached_to_no_iceberg = (mask.e == ICEBERGMASK_NO_ICEBERG ||
                                    mask.w == ICEBERGMASK_NO_ICEBERG ||
                                    mask.n == ICEBERGMASK_NO_ICEBERG ||
                                    mask.s == ICEBERGMASK_NO_ICEBERG);

        if (vIcebergMask(i, j) == ICEBERGMASK_ICEBERG_CAND &&
            (attached_to_grounded || attached_to_no_iceberg)) {

          vIcebergMask(i, j) = ICEBERGMASK_NO_ICEBERG;
          done = false;
        }

      }
    }
    ierr = vIcebergMask.end_access(); CHKERRQ(ierr);

    ierr = vIcebergMask.update_ghosts(); CHKERRQ(ierr);

    // We're "done" only if the iceberg mask stopped changing on *all*
    // processor sub-domains.
    int flag = done;
    MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_LAND, grid.com);
    done = flag;

    loopcount += 1;
  }

  ierr = verbPrintf(3, grid.com,
    "PISM-PIK INFO:  %d loop(s) were needed to identify whether there are icebergs \n",
    loopcount); CHKERRQ(ierr);

  return 0;
}


/*!
 * We have distinguished icebergs from attached floating regions in
 * identifyNotAnIceBerg(). Now we actually eliminate the former (meaning we set
 * the ice thickness to zero and mark the boxes as icefree-ocean) and leave the
 * latter as they are.
 */
PetscErrorCode IceModel::killIdentifiedIceBergs() {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "######### killIdentifiedIceBergs() start\n"); CHKERRQ(ierr);
  const bool vpik = config.get_flag("verbose_pik_messages");
  
  PetscScalar
    my_discharge_flux = 0,
    discharge_flux = 0; 
  const PetscScalar dx = grid.dx, dy = grid.dy;

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vh.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      if (vIcebergMask(i, j) == ICEBERGMASK_ICEBERG_CAND) {
         // it is not a candidate any more, it is an iceberg!
         my_discharge_flux -= vH(i, j);
         vH(i, j) = 0.0;
         vh(i, j) = 0.0;
         vMask(i, j) = MASK_ICE_FREE_OCEAN;
         if (vpik) {
           PetscSynchronizedPrintf(grid.com, 
                     "PISM-PIK INFO: [rank %d] killed iceberg at i=%d, j=%d\n", 
                     grid.rank, i, j);
         }
      }

    }
  }

  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);
  
  ierr = PISMGlobalSum(&my_discharge_flux,     &discharge_flux,     grid.com); CHKERRQ(ierr);
  PetscScalar factor = config.get("ice_density") * (dx * dy);
  discharge_flux_cumulative     += discharge_flux     * factor;

  ierr = vMask.update_ghosts(); CHKERRQ(ierr);
  ierr = vH.update_ghosts(); CHKERRQ(ierr);
  ierr = vh.update_ghosts(); CHKERRQ(ierr);

  return 0;
}


/*!
 * This routine is used when strain-rate based calving is applied. It kills
 * single grid cell icebergs and one-grid cell wide ice 'noses', because no
 * proper strain rate eigenvalues can be derived there.
 */
PetscErrorCode IceModel::killEasyIceBergs() {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "######### killEasyIceBergs() start\n"); CHKERRQ(ierr);
  const bool vpik = config.get_flag("verbose_pik_messages");
  
  PetscScalar
    my_discharge_flux = 0,
    discharge_flux = 0;
  const PetscScalar dx = grid.dx, dy = grid.dy;

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);

  PetscReal sea_level;
  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else { SETERRQ(grid.com, 2, "PISM ERROR: ocean == NULL"); }

  double ocean_rho = config.get("sea_water_density"),
    ice_rho = config.get("ice_density");

  // looking for grid-cell wide floating ice noses that have at least six neighbors
  // of thickness H=0 like this (o ocean, fl and x floating):
  //
  // o  o  o        o o o        o  o  o
  // fl x fl   OR   o x o   OR   o  x fl
  // o  o  o        o o o        o  o  o

  PetscReal C = (1.0 - ice_rho / ocean_rho);

  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vh.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      planeStar<PetscScalar> thk = vH.star(i, j),
        bed = vbed.star(i, j);

      // instead of updating surface elevation, counting here floating or icefree neighbors
      const PetscScalar hgrounded = bed.ij + thk.ij,
        hfloating = sea_level + (1.0 - ice_rho / ocean_rho) * thk.ij;

      if (vH(i, j) > 0.0 && hgrounded < hfloating) { //is floating ice shelf

        const PetscScalar
          hgrounded_eb = bed.e + thk.e,
          hfloating_eb = sea_level + C * thk.e,
          hgrounded_wb = bed.w + thk.w,
          hfloating_wb = sea_level + C * thk.w,
          hgrounded_nb = bed.n + thk.n,
          hfloating_nb = sea_level + C * thk.n,
          hgrounded_sb = bed.s + thk.s,
          hfloating_sb = sea_level + C * thk.s;

        PetscInt jcount = 0, icount = 0; // grid-cell wide floating ice nose

        if (vH(i + 1, j + 1) == 0.0) {jcount += 1; icount += 1;}
        if (vH(i + 1, j - 1) == 0.0) {jcount += 1; icount += 1;}
        if (vH(i - 1, j + 1) == 0.0) {jcount += 1; icount += 1;}
        if (vH(i - 1, j - 1) == 0.0) {jcount += 1; icount += 1;}

        if (thk.e == 0.0) jcount += 1;
        if (thk.w == 0.0) jcount += 1;
        if (thk.n == 0.0) icount += 1;
        if (thk.s == 0.0) icount += 1;

        if ((icount == 6 && hgrounded_eb < hfloating_eb && hgrounded_wb < hfloating_wb) ||
            (jcount == 6 && hgrounded_nb < hfloating_nb && hgrounded_sb < hfloating_sb)) {

	  my_discharge_flux -= vHnew(i, j);
          vHnew(i, j) = 0.0;
          vh(i, j) = 0.0;
          vMask(i, j) = MASK_ICE_FREE_OCEAN;
          if (vpik) {
            PetscSynchronizedPrintf(grid.com, 
              "PISM-PIK INFO: [rank %d] cut off nose or one-box-iceberg at i=%d, j=%d\n",
              grid.rank, i, j);
          }
        }
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);

  // finally copy vHnew into vH and communicate ghosted values
  ierr = vHnew.update_ghosts(vH); CHKERRQ(ierr);

  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);

  // looking for one-grid-cell icebergs, that have 4 neighbors of thickness H=0
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      planeStar<PetscScalar> thk = vH.star(i, j);

      // instead of updating surface elevation, counting here floating or icefree neighbors
      const PetscScalar hgrounded = vbed(i, j) + thk.ij,
        hfloating = sea_level + C * thk.ij;

      bool all_4neighbors_icefree = (thk.e == 0.0 && thk.w == 0.0 &&
                                     thk.n == 0.0 && thk.s == 0.0);

      if (thk.ij > 0.0 && hgrounded < hfloating && all_4neighbors_icefree) {
	my_discharge_flux -= vHnew(i, j);
        vHnew(i, j) = 0.0;
        vh(i, j) = 0.0;
        vMask(i, j) = MASK_ICE_FREE_OCEAN;
        if (vpik) {
          PetscSynchronizedPrintf(grid.com,
            "PISM-PIK INFO: [rank %d] killed isolated one-box-iceberg at i=%d, j=%d\n",
            grid.rank, i, j);
        }
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);

  ierr = vHnew.update_ghosts(vH); CHKERRQ(ierr);

  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);

  // looking for one-grid-cell partially filled grid cells, that have 4 neighbors of thickness H=0
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vHref.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      // instead of updating surface elevation, counting here floating or icefree neighbors
      bool all_4neighbors_icefree = (vH(i + 1, j) == 0.0 &&
                                     vH(i - 1, j) == 0.0 &&
                                     vH(i, j + 1) == 0.0 &&
                                     vH(i, j - 1) == 0.0);
      if (vHref(i, j) > 0.0 && all_4neighbors_icefree) {
	my_discharge_flux -= vHref(i, j);
        vHref(i, j) = 0.0;
        if (vpik) {
          PetscSynchronizedPrintf(grid.com, 
            "PISM-PIK INFO: [rank %d] killed lonely partially filled grid cell at i = %d, j = %d\n",
            grid.rank, i, j);
        }
      }
    }
  }
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);
  ierr = vHref.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  
  ierr = PISMGlobalSum(&my_discharge_flux,     &discharge_flux,     grid.com); CHKERRQ(ierr);
  PetscScalar factor = config.get("ice_density") * (dx * dy);
  discharge_flux_cumulative     += discharge_flux     * factor;

  if (vpik)  // actually get output from PetscSynchronizedPrintf()
    PetscSynchronizedFlush(grid.com);

  ierr = vHnew.update_ghosts(vH); CHKERRQ(ierr);
  ierr = vMask.update_ghosts(); CHKERRQ(ierr);
  ierr = vh.update_ghosts(); CHKERRQ(ierr);
  ierr = vHref.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

