// Copyright (C) 2011, 2012 PISM Authors
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

#include "iceModel.hh"
#include "IceGrid.hh"

//! \file iMhydrology.cc Currently, only the most minimal possible hydrology model: diffusion of stored basal water.

//! Explicit time step for diffusion of subglacial water layer bwat.
/*!
See equation (11) in \ref BBssasliding , namely
  \f[W_t = K \nabla^2 W.\f]
The diffusion constant \f$K\f$ is chosen so that the fundamental solution (Green's
function) of this equation has standard deviation \f$\sigma=L\f$ at time t=\c diffusion_time.
Note that \f$2 \sigma^2 = 4 K t\f$.

The time step restriction for the explicit method for this equation is believed
to be so rare, for most values of \c bwat_diffusion_distance and \c bwat_diffusion_time
that we simply halt execution if it occurs.

Uses vWork2d[0] to temporarily store new values for bwat.
 */
PetscErrorCode IceModel::diffuse_bwat() {
  PetscErrorCode  ierr;

  const PetscScalar
    L = config.get("bwat_diffusion_distance"),
    diffusion_time = config.get("bwat_diffusion_time", "years", "seconds"); // convert to seconds

  const PetscScalar K = L * L / (2.0 * diffusion_time),
                    Rx = K * dt_TempAge / (grid.dx * grid.dx),
                    Ry = K * dt_TempAge / (grid.dy * grid.dy),
                    oneM4R = 1.0 - 2.0 * Rx - 2.0 * Ry;
  if (oneM4R <= 0.0) {
    SETERRQ(grid.com, 1,
       "PISM ERROR: diffuse_bwat() requires 1 - 2Rx - 2Ry <= 0 and an explicit step;\n"
       "  current timestep would cause this method to be unstable; this event is believed\n"
       "  to be so rare that the timestep restriction is not part of the adaptive scheme\n"
       "ENDING!\n\n");
  }

  // communicate ghosted values so neighbors are valid;
  // note that temperatureStep() and enthalpyAndDrainageStep() modify vbwat,
  // but they do not update ghosts because only the current process needs that
  ierr = vbwat.beginGhostComm(); CHKERRQ(ierr);
  ierr = vbwat.endGhostComm(); CHKERRQ(ierr);

  PetscScalar **bwatnew; 
  ierr = vbwat.begin_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(bwatnew); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      bwatnew[i][j] = oneM4R * vbwat(i,j)
                       + Rx * (vbwat(i+1,j  ) + vbwat(i-1,j  ))
                       + Ry * (vbwat(i  ,j+1) + vbwat(i  ,j-1));
    }
  }
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vbwat.end_access(); CHKERRQ(ierr);

  // finally copy new into vbwat and communicate ghosts at the same time
  ierr = vWork2d[0].beginGhostComm(vbwat); CHKERRQ(ierr);
  ierr = vWork2d[0].endGhostComm(vbwat); CHKERRQ(ierr);

  return 0;
}

