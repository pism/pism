// Copyright (C) 2007--2008 Jed Brown and Ed Bueler
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

// this code fragment is for inclusion in IceModel::dumpToFile_netCDF()
// it was *not* automatically generated

  if (grid.rank == 0) {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "history",
                           strlen(history), history);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
    
    stat = nc_enddef (ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    double t = grid.year * secpera;
    size_t zero = 0;
    stat = nc_put_var1_double(ncid, t_id, &zero, &t); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    ierr = nct.put_dimension_regular(ncid, x_id, x_len, -grid.Lx, grid.dx); CHKERRQ(ierr);
    ierr = nct.put_dimension_regular(ncid, y_id, y_len, -grid.Ly, grid.dy); CHKERRQ(ierr);
    
    ierr = nct.put_dimension(ncid, z_id, z_len, grid.zlevels); CHKERRQ(ierr);
    ierr = nct.put_dimension(ncid, zb_id, zb_len, grid.zblevels); CHKERRQ(ierr);

    stat = nc_put_att_double(ncid, polar_stereographic_id, "straight_vertical_longitude_from_pole",
                            NC_DOUBLE, 1, &psParams.svlfp); CHKERRQ(nc_check(stat));
    stat = nc_put_att_double(ncid, polar_stereographic_id, "latitude_of_projection_origin",
                            NC_DOUBLE, 1, &psParams.lopo); CHKERRQ(nc_check(stat));
    stat = nc_put_att_double(ncid, polar_stereographic_id, "standard_parallel",
                            NC_DOUBLE, 1, &psParams.sp); CHKERRQ(nc_check(stat));

    if (useSSAVelocity) {
      int one = 1;
      stat = nc_put_att_int(ncid, NC_GLOBAL, "haveSSAvelocities",
			    NC_INT, 1, &one); CHKERRQ(nc_check(stat));
    }
  }

  // 2-D model quantities
  ierr = nct.put_local_var(ncid, lon_id, grid.da2, vLongitude, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, lat_id, grid.da2, vLatitude, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, mask_id, grid.da2, vMask, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, thk_id, grid.da2, vH, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, bwat_id, grid.da2, vHmelt, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, topg_id, grid.da2, vbed, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, dbdt_id, grid.da2, vuplift, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  if (useSSAVelocity) {
    ierr = nct.put_local_var(ncid, vubarSSA_id, grid.da2, vubarSSA, g2,
			     s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
    ierr = nct.put_local_var(ncid, vvbarSSA_id, grid.da2, vvbarSSA, g2,
			     s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  }

  // 3-D model quantities
  ierr = T3.putVecNC(ncid, s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = Tb3.putVecNC(ncid, s, cb, 4, a_mpi, max_a_len); CHKERRQ(ierr);
//  ierr = nct.put_local_var(&grid, ncid, litho_temp_id, grid.da3b, Tb3.v, g3b,
//                       s, cb, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = tau3.putVecNC(ncid, s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  // 2-D climate quantities
  ierr = nct.put_local_var(ncid, artm_id, grid.da2, vTs, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, bheatflx_id, grid.da2, vGhf, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, acab_id, grid.da2, vAccum, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
      // write tillphi = till friction angle in degrees
  ierr = nct.put_local_var(ncid, tillphi_id, grid.da2, vtillphi, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // 2-D diagnostic quantities
  // note h is diagnostic because it is recomputed by h=H+b at each time step
  // these are not written in MKS units because they are intended to be viewed,
  // not read by programs; IS THIS THE RIGHT CHOICE?
  ierr = nct.put_local_var(ncid, usurf_id, grid.da2, vh, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = VecCopy(vdHdt,vWork2d[0]); CHKERRQ(ierr);
  ierr = VecScale(vWork2d[0],secpera); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, dHdt_id, grid.da2, vWork2d[0], g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // compute cbar = sqrt(ubar^2 + vbar^2) and save it
  ierr = VecPointwiseMult(vWork2d[0], vubar, vubar); CHKERRQ(ierr);
  ierr = VecPointwiseMult(vWork2d[1], vvbar, vvbar); CHKERRQ(ierr);
  ierr = VecAXPY(vWork2d[0], 1, vWork2d[1]); CHKERRQ(ierr);
  PetscScalar **a, **H;
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &a); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        a[i][j] = sqrt(a[i][j]) * secpera;
      } else {
        a[i][j] = 0.0; // no ice at location
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &a); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, cbar_id, grid.da2, vWork2d[0], g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // compute cflx = cbar .* thk and save it
  ierr = VecPointwiseMult(vWork2d[1], vWork2d[0], vH); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, cflx_id, grid.da2, vWork2d[1], g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // compute cbase  = sqrt(u|_{z=0}^2 + v|_{z=0}^2) and save it
  ierr = u3.needAccessToVals(); CHKERRQ(ierr);
  ierr = v3.needAccessToVals(); CHKERRQ(ierr);
  ierr = u3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);
  ierr = v3.getHorSlice(vWork2d[1], 0.0); CHKERRQ(ierr);
  ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
  PetscScalar **ub, **vb;
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &a); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      if (H[i][j] > 0.0) {
        a[i][j] = sqrt(PetscSqr(ub[i][j]) + PetscSqr(vb[i][j])) * secpera;
      } else {
        a[i][j] = 0.0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &a); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, cbase_id, grid.da2, vWork2d[2], g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // compute csurf = sqrt(u|_surface^2 + v|_surface^2) and save it
  ierr = u3.needAccessToVals(); CHKERRQ(ierr);
  ierr = v3.needAccessToVals(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValuesVec2d(vWork2d[0], vH); CHKERRQ(ierr);
  ierr = v3.getSurfaceValuesVec2d(vWork2d[1], vH); CHKERRQ(ierr);
  ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
  PetscScalar **us, **vs;
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &us); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &vs); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &a); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      if (H[i][j] > 0.0) {
        a[i][j] = sqrt(PetscSqr(us[i][j]) + PetscSqr(vs[i][j])) * secpera;
      } else {
        a[i][j] = 0.0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &us); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &vs); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &a); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, csurf_id, grid.da2, vWork2d[2], g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // compute wsurf, the surface values of vertical velocity
  ierr = w3.needAccessToVals(); CHKERRQ(ierr);
  ierr = w3.getSurfaceValuesVec2d(vWork2d[0], vH); CHKERRQ(ierr);
  ierr = w3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = VecScale(vWork2d[0],secpera); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, wsurf_id, grid.da2, vWork2d[0], g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // compute magnitude of basal shear stress = rho g H |grad h|
  ierr = computeDrivingStress(vWork2d[0],vWork2d[1]); CHKERRQ(ierr);
  ierr = getMagnitudeOf2dVectorField(vWork2d[0],vWork2d[1],vWork2d[2]); CHKERRQ(ierr);
  ierr = nct.put_local_var(ncid, taud_id, grid.da2, vWork2d[2], g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // write out yield stress
  ierr = nct.put_local_var(ncid, tauc_id, grid.da2, vtauc, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

