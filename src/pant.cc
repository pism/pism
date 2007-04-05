// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

static char help[] =
  "PISM driver whose purpose is to have either: (1) spatially-variable scalar basal sliding\n"
  "friction parameter beta, (2) spatially-variable vector basal sliding friction parameter\n"
  "(beta_x,beta_y), or (3) a Schoof-type plastic till model with dependence on a\n"
  "very basic basal water pressure model.  The spatially-varying beta (or vector version)\n"
  "will derive from balance velocities.\n";

/* 

1. For distributed *vector* beta computed from balance velocities:

     pant -if ant40km.nc -mv -gk -e 1.2 -y 10 -verbose -balvel init.nc -betaxy betaxymap

Here ant40km.nc contains a saved model state with a vaguely 
reasonable temperature field.  The balance velocities are read in from init.nc.
Then the deformation (SIA) horizontal velocity is computed and subtracted from
balance velocities to get a sliding velocity field (every where on grounded ice).
The MacAyeal equations are (i.e. the MacAyeal nonlinear operator is) applied to 
determine the vector (beta_x,beta_y).  These are written to a NetCDF file betaxymap.nc.
Then a 10 year run is done using these; the model should be nearly in steady state 
(although it will drift away as the temperature changes).

2. For distributed *scalar* beta computed from balance velocities:

     pant -if ant40km.nc -mv -gk -e 1.2 -y 10 -verbose -balvel init.nc -beta betamap

[CONCRETELY:
     pant -if ant150k_40kmSMMV.nc -mv -gk -e 1.2 -ksp_rtol 1e-6 -ocean_kill -verbose\   
         -balvel init.nc -beta betamap -d Bncm -y 10
]

As in 1. except that once the vector (beta_x,beta_y) is generated we combine: beta
= (u^2 beta_x + v^2 beta_y)/(u^2+v^2) with max of BETAMX = 1.0e13.  beta is written to 
a NetCDF file betamap.nc.  Then a 10 year run is done using the distributed beta; 
the model should less in steady state than in 1.

3. For Schoof-type solution of plastic till problem, with yield stress computed from 
thermal model and estimate of basal till pore water pressure:

     pant -if ant40km.nc -mv -gk -e 1.2 -y 10 -verbose -plastic taucmap

Here the plastic till version of the MacAyeal equations are solved from the saved model
state ant40km.nc.  In particular, the vHmelt values and the basal temps (vT) from
ant40km.nc are used to comput vtauc.  This map is saved in taucmap.nc.  Then 10 years
are run.

*/

#include <cstring>
#include <petscbag.h>

#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

#include <netcdf.h>
#include "nc_util.hh"

class IceDragYieldModel : public IceModel {
public:
  PetscTruth     balvelFileSet, doBetaxy, doPlastic;
  char           balvelFileName[PETSC_MAX_PATH_LEN], betaFileName[PETSC_MAX_PATH_LEN],
                 betaxyFileName[PETSC_MAX_PATH_LEN], taucFileName[PETSC_MAX_PATH_LEN];
  Vec            vbetax, vbetay;
  Vec            vmagbalvel, vNu[2];
  Vec 	         dragxy;  // this is not da2 sized; it's a (2 Mx My) x 1 column
  
  IceDragYieldModel(IceGrid &g, IceType &i);
  PetscErrorCode createDragYieldVecs();
  PetscErrorCode destroyDragYieldVecs();
  PetscErrorCode dragYieldOptions();
  
  // for distribute beta (and (beta_x,beta_y) both)
  PetscErrorCode computeScalarBeta();
  PetscErrorCode readBalVels();
  PetscErrorCode getEffectiveViscosity();
  PetscErrorCode computeDragFromBalanceVelocity();
  PetscErrorCode moveDragxytoDAbetas();
  PetscErrorCode writeBetaNC();
  PetscErrorCode writeMatlabDragFile(const char *basename);
  PetscErrorCode updateMaskFromBeta();
//  PetscErrorCode computeVectorBeta();
//  PetscErrorCode writeBetaxyNC();
  
  // for plastic case
//  PetscErrorCode computePlasticTauc();
//  PetscErrorCode writeTaucNC();

};


#define BETAMAX 1.0e13  // Pa s m^-1; max allowed drag coeff beta


IceDragYieldModel::IceDragYieldModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {
  // do nothing (except pass g,i to constructor for IceModel)
}


PetscErrorCode IceDragYieldModel::createDragYieldVecs() {
  PetscErrorCode ierr;
  
  ierr = VecDuplicate(vh, &vbetax); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vbetay); CHKERRQ(ierr);
  const PetscInt M = 2 * grid.p->Mx * grid.p->My;
  ierr = VecCreateMPI(grid.com, PETSC_DECIDE, M, &dragxy); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragYieldModel::destroyDragYieldVecs() {
  PetscErrorCode  ierr;

  ierr = VecDestroy(vbetax); CHKERRQ(ierr);
  ierr = VecDestroy(vbetay); CHKERRQ(ierr);
  ierr = VecDestroy(dragxy); CHKERRQ(ierr);
  // note vmagbalvel, vNu[0], vNu[1] are just IceModel work Vecs
  return 0;
}


PetscErrorCode IceDragYieldModel::dragYieldOptions() {
  PetscErrorCode  ierr;
  PetscTruth      plasticSet, betaSet, betaxySet;

  // get filenames
  ierr = PetscOptionsGetString(PETSC_NULL, "-balvel", balvelFileName,
                               PETSC_MAX_PATH_LEN, &balvelFileSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-beta", betaFileName,
                               PETSC_MAX_PATH_LEN, &betaSet); CHKERRQ(ierr);
  strcat(betaFileName, ".nc");

  ierr = PetscOptionsHasName(PETSC_NULL, "-plastic", &plasticSet); CHKERRQ(ierr);  // filename ignored
  ierr = PetscOptionsHasName(PETSC_NULL, "-betaxy", &betaxySet); CHKERRQ(ierr);  // filename ignored

  if (plasticSet == PETSC_TRUE) {
    if (betaSet == PETSC_TRUE) {
      SETERRQ(1,"option conflict in pant; both -plastic and -beta set");
    }
    if (betaxySet == PETSC_TRUE) {
      SETERRQ(2,"option conflict in pant; both -plastic and -betaxy set");
    }
    if (balvelFileSet == PETSC_TRUE) {
      SETERRQ(3,"option conflict in pant; both -plastic and a balance velocities file (-balvel) set");
    }
    doPlastic = PETSC_TRUE;
    doBetaxy = PETSC_FALSE;
  } else {
    if (balvelFileSet == PETSC_FALSE) {
      SETERRQ(4,"pant requires source for balance velocities if -plastic is not set; use -balvel foo.nc");
    }
    doPlastic = PETSC_FALSE;
    if (betaSet == PETSC_TRUE) {
      if (betaxySet == PETSC_TRUE) {
        SETERRQ(5,"option conflict in pant; both -beta and -betaxy set");
      }
      doBetaxy = PETSC_FALSE;
    } else {
      doBetaxy = PETSC_TRUE;
    }
  }
  return 0;
}


PetscErrorCode nc_check(int stat) {
  if (stat)
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
  return 0;
}


PetscErrorCode IceDragYieldModel::readBalVels() {
  PetscErrorCode  ierr;
  int stat, ncid, v_balvel;

  if (grid.rank == 0) {
    stat = nc_open(balvelFileName, 0, &ncid); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "balvel", &v_balvel); CHKERRQ(nc_check(stat));
  }

  vmagbalvel = vWork2d[2];  // already allocated
  Vec vzero;
  VecScatter ctx;
  ierr = VecScatterCreateToZero(g2, &ctx, &vzero); CHKERRQ(ierr);  
  ierr = getIndZero(grid.da2, g2, vzero, ctx); CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_balvel, grid.da2, vmagbalvel, g2, vzero); CHKERRQ(ierr);
  ierr = VecDestroy(vzero); CHKERRQ(ierr);
  ierr = VecScatterDestroy(ctx); CHKERRQ(ierr);
  ierr = VecScale(vmagbalvel,1.0/secpera); CHKERRQ(ierr);  // convert to m/s

  if (grid.rank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }
  return 0;
}


PetscErrorCode IceDragYieldModel::writeBetaNC() {
  PetscErrorCode  ierr;

/* BLOCK of code here parallels write_attributes.c */

   int  stat;			/* return status */
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int x_dim;
   int y_dim;
   int t_dim;

   /* dimension lengths */
   size_t x_len = grid.p->Mx;
   size_t y_len = grid.p->My;
   size_t t_len = NC_UNLIMITED;

   /* variable ids */
   int x_id;
   int y_id;
   int t_id;
   int beta_id;
   int betax_id;
   int betay_id;

#  define RANK_x 1
#  define RANK_y 1
#  define RANK_t 1
#  define RANK_beta 3
#  define RANK_betax 3
#  define RANK_betay 3

   /* variable shapes */
   int x_dims[RANK_x];
   int y_dims[RANK_y];
   int t_dims[RANK_t];
   int beta_dims[RANK_beta];
   int betax_dims[RANK_betax];
   int betay_dims[RANK_betay];

   /* enter define mode */
if (grid.rank == 0) {
   stat = nc_create(betaFileName, NC_CLOBBER|NC_64BIT_OFFSET, &ncid);
   CHKERRQ(nc_check(stat));

   /* define dimensions */
   stat = nc_def_dim(ncid, "x", x_len, &x_dim);
   CHKERRQ(nc_check(stat));
   stat = nc_def_dim(ncid, "y", y_len, &y_dim);
   CHKERRQ(nc_check(stat));
   stat = nc_def_dim(ncid, "t", t_len, &t_dim);
   check_err(stat,__LINE__,__FILE__);

   /* define variables */

   x_dims[0] = x_dim;
   stat = nc_def_var(ncid, "x", NC_FLOAT, RANK_x, x_dims, &x_id);
   CHKERRQ(nc_check(stat));

   y_dims[0] = y_dim;
   stat = nc_def_var(ncid, "y", NC_FLOAT, RANK_y, y_dims, &y_id);
   CHKERRQ(nc_check(stat));

   t_dims[0] = t_dim;
   stat = nc_def_var(ncid, "t", NC_FLOAT, RANK_t, t_dims, &t_id);
   check_err(stat,__LINE__,__FILE__);

   beta_dims[0] = t_dim;
   beta_dims[1] = x_dim;
   beta_dims[2] = y_dim;
   stat = nc_def_var(ncid, "beta", NC_FLOAT, RANK_beta, beta_dims, &beta_id);
   CHKERRQ(nc_check(stat));

   betax_dims[0] = t_dim;
   betax_dims[1] = x_dim;
   betax_dims[2] = y_dim;
   stat = nc_def_var(ncid, "betax", NC_FLOAT, RANK_betax, betax_dims, &betax_id);
   CHKERRQ(nc_check(stat));

   betay_dims[0] = t_dim;
   betay_dims[1] = x_dim;
   betay_dims[2] = y_dim;
   stat = nc_def_var(ncid, "betay", NC_FLOAT, RANK_betay, betay_dims, &betay_id);
   CHKERRQ(nc_check(stat));

   /* assign attributes */
   stat = nc_put_att_text(ncid, x_id, "axis", 1, "X");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, x_id, "long_name", 32, "x-coordinate in Cartesian system");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, x_id, "standard_name", 23, "projection_x_coordinate");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, x_id, "units", 1, "m");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, y_id, "axis", 1, "Y");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, y_id, "long_name", 32, "y-coordinate in Cartesian system");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, y_id, "standard_name", 23, "projection_y_coordinate");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, y_id, "units", 1, "m");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, t_id, "long_name", 4, "time");
   check_err(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, t_id, "units", 33, "seconds since 2007-01-01 00:00:00");
   check_err(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, t_id, "calendar", 4, "none");
   check_err(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, t_id, "axis", 1, "T");
   check_err(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, beta_id, "long_name", 30, "ice_basal_friction_coefficient");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, beta_id, "units", 8, "Pa s m-1");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, betax_id, "long_name", 42, "ice_basal_friction_coefficient_x_component");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, betax_id, "units", 8, "Pa s m-1");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, betay_id, "long_name", 42, "ice_basal_friction_coefficient_y_component");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, betay_id, "units", 8, "Pa s m-1");
   CHKERRQ(nc_check(stat));

} // end if (grid.rank == 0)
/* ENDBLOCK */


  if (grid.rank == 0) {
    stat = nc_enddef(ncid);
    CHKERRQ(nc_check(stat));

    float t = grid.p->year * secpera;
    size_t zero = 0;
    stat = nc_put_var1_float(ncid, t_id, &zero, &t);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    ierr = put_dimension_regular(ncid, x_id, x_len, -grid.p->Lx, grid.p->dx); CHKERRQ(ierr);
    ierr = put_dimension_regular(ncid, y_id, y_len, -grid.p->Ly, grid.p->dy); CHKERRQ(ierr);
  }

  int s[] = {0, grid.xs, grid.ys, 0};            // Start local block: t dependent
  int c[] = {1, grid.xm, grid.ym, grid.p->Mz};   // Count local block: t dependent
 
  // Allocate some memory. 
  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym;
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(float), &a_mpi); CHKERRQ(ierr);

  // 2-D model quantities
  ierr = put_local_var(&grid, ncid, beta_id, NC_FLOAT, grid.da2, vbeta, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, betax_id, NC_FLOAT, grid.da2, vbetax, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, betay_id, NC_FLOAT, grid.da2, vbetay, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // We are done with these buffers
  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_close(ncid);
    CHKERRQ(nc_check(stat));
  }

  return 0;
}


PetscErrorCode IceDragYieldModel::getEffectiveViscosity() {
  PetscErrorCode ierr;
  PetscScalar epsilon = 0.0;  // actually get effective viscosity (not bdd below)
  PetscScalar **mask;
  
  // need effective viscosity *everywhere* so mark all grounded points as DRAGGING
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) != MASK_FLOATING)   mask[i][j] = MASK_DRAGGING;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  // communicate, because computing eff visc requires neighbor vals of mask
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  // put effective viscosity in working arrays
  vNu[0] = vWork2d[0];
  vNu[1] = vWork2d[1];
  ierr = computeEffectiveViscosity(vNu, epsilon); CHKERRQ(ierr);  

  // view if  "-d n"
  if (nuView[0] != PETSC_NULL && nuView[1] != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vNu[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, nuView[0]); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vNu[1], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, nuView[1]); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com,"\n final viscosity; pausing for 10 seconds ... \n"); CHKERRQ(ierr);
    ierr = PetscSleep(10); CHKERRQ(ierr);
  }
  return 0;
}

  
// compute drag from velocities and mass-balance-velocity
PetscErrorCode IceDragYieldModel::computeDragFromBalanceVelocity() {
  PetscErrorCode ierr;
  Mat A;
  Vec x, result, rhs;
  PetscScalar **u, **v, **mask, **balvel, **beta, **tauc;

  // allocate Vecs
  const PetscInt M = 2 * grid.p->Mx * grid.p->My;
  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE, M, M,
                         13, PETSC_NULL, 13, PETSC_NULL,
                         &A); CHKERRQ(ierr);
  ierr = VecDuplicate(dragxy, &result); CHKERRQ(ierr);
  ierr = VecDuplicate(dragxy, &rhs); CHKERRQ(ierr);
  ierr = VecDuplicate(dragxy, &x); CHKERRQ(ierr);

  // build discrete version of MacAyeal-Morland equations (A x = rhs) at all grounded points
  ierr = assembleMacayealMatrix(vNu, A); CHKERRQ(ierr);
  ierr = assembleMacayealRhs(false, rhs); CHKERRQ(ierr);

  // Note rhs contains driving terms  \rho g H \grad h  but  A  contains basal
  // drag term [betax*u betay*v]'.  It must be removed.
  // Also set x = [u v]'  (i.e. interleaved).
  ierr = VecSet(x, 0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vmagbalvel, &balvel); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt    J = 2*j;
      const PetscInt    rowU = i*2*grid.p->My + J;
      const PetscInt    rowV = i*2*grid.p->My + J+1;
      // remove old drag term
      if (intMask(mask[i][j]) == MASK_DRAGGING) {
        ierr = MatSetValue(A, rowU, rowU,
                           - basalDragx(beta, tauc, u, v, i, j), ADD_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(A, rowV, rowV,
                           - basalDragy(beta, tauc, u, v, i, j), ADD_VALUES); CHKERRQ(ierr);
      }
      // remove deformational from mass-balance velocities; put in x
      const PetscScalar c = sqrt(PetscSqr(u[i][j]) + PetscSqr(v[i][j]));
      PetscScalar       residualFactor = (balvel[i][j] / c) - 1.0;
      // if (residualFactor < 0.0)   residualFactor = 0.0;  // avoid negative drag coeff here?  probably not
      ierr = VecSetValue(x, rowU, residualFactor * u[i][j], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(x, rowV, residualFactor * v[i][j], INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vmagbalvel, &balvel); CHKERRQ(ierr);
  // communicate!
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

  // With new A and x, compute  result = Ax - rhs.
  ierr = MatMult(A,x,result); CHKERRQ(ierr);
  ierr = VecAXPY(result,-1.0,rhs); CHKERRQ(ierr);  // note result = 0 if MASK_SHEET
  
  // As result is [betax*u betay*v]', divide by velocities, but only where grounded!
  // where floating, report no drag coeff
  ierr = VecPointwiseDivide(dragxy,result,x); CHKERRQ(ierr);
  ierr = VecScale(dragxy,-1.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt    J = 2*j;
      const PetscInt    rowU = i*2*grid.p->My + J;
      const PetscInt    rowV = i*2*grid.p->My + J+1;
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        ierr = VecSetValue(dragxy, rowU, 0.0, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(dragxy, rowV, 0.0, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dragxy); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dragxy); CHKERRQ(ierr);
  
  ierr = MatDestroy(A); CHKERRQ(ierr);
  ierr = VecDestroy(x); CHKERRQ(ierr);
  ierr = VecDestroy(result); CHKERRQ(ierr);
  ierr = VecDestroy(rhs); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragYieldModel::moveDragxytoDAbetas() {
  PetscErrorCode  ierr;
  PetscScalar     **u, **v, **betax, **betay, **beta, *betaxylong;
  Vec             dragxyLoc;
  VecScatter      dragxyScatterGlobalToLocal;

  // Move the solution onto a grid which can be accessed normally. Since the parallel
  // layout of the vector dragxy does not in general have anything to do with the DA based
  // vectors, we must scatter the entire vector to all processors.
  const PetscInt MM = 2 * grid.p->Mx * grid.p->My;
  ierr = VecCreateSeq(PETSC_COMM_SELF, MM, &dragxyLoc);
//  ierr = VecScatterCreate(dragxy, PETSC_NULL, dragxyLoc, PETSC_NULL,
//                          &dragxyScatterGlobalToLocal); CHKERRQ(ierr);
  ierr = VecScatterCreateToAll(dragxy, &dragxyScatterGlobalToLocal, &dragxyLoc); CHKERRQ(ierr);
  ierr = VecScatterBegin(dragxy, dragxyLoc, INSERT_VALUES, SCATTER_FORWARD,
                         dragxyScatterGlobalToLocal); CHKERRQ(ierr);
  ierr = VecScatterEnd(dragxy, dragxyLoc, INSERT_VALUES, SCATTER_FORWARD,
                       dragxyScatterGlobalToLocal); CHKERRQ(ierr);

  /* now extract betax,betay and then compute average beta */
  ierr = VecGetArray(dragxyLoc, &betaxylong); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbetax, &betax); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbetay, &betay); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  const PetscInt     M = 2 * grid.p->My;
  const PetscScalar  BETAMIN = 1.0e6;
  const PetscScalar  epsuv = PetscSqr(1.0 / secpera);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
    
      betax[i][j] = betaxylong[i*M + 2*j];
      if (isnan(betax[i][j]))  betax[i][j] = BETAMAX;
      if (isinf(PetscAbs(betax[i][j])))  betax[i][j] = BETAMAX;
      betay[i][j] = betaxylong[i*M + 2*j+1];
      if (isnan(betay[i][j]))  betay[i][j] = BETAMAX;
      if (isinf(PetscAbs(betay[i][j])))  betay[i][j] = BETAMAX;

      const PetscScalar usqr = PetscSqr(u[i][j])+epsuv, vsqr = PetscSqr(v[i][j])+epsuv; 
      beta[i][j] = (usqr * betax[i][j] + vsqr * betay[i][j]) / (usqr + vsqr);
      if (beta[i][j] > BETAMAX)  beta[i][j] = BETAMAX;
      if (beta[i][j] < BETAMIN)  beta[i][j] = BETAMAX;  // this is better value for future use
    }
  }
  ierr = VecRestoreArray(dragxyLoc, &betaxylong); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbetax, &betax); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbetay, &betay); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);

  ierr = VecDestroy(dragxyLoc); CHKERRQ(ierr);
  ierr = VecScatterDestroy(dragxyScatterGlobalToLocal); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragYieldModel::updateMaskFromBeta() {
  PetscErrorCode  ierr;
  PetscScalar     **mask, **beta;

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) != MASK_FLOATING) {
        if (beta[i][j] < 1e-2 * BETAMAX) {
          mask[i][j] = MASK_DRAGGING;
        } else {
          mask[i][j] = MASK_SHEET;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  // now remove singleton and nearly singleton DRAGGING points, i.e. if 
  //   all BOX stencil neighbors are SHEET or if at most one is DRAGGING:
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (intMask(mask[i][j]) == MASK_DRAGGING) {
        const int neighmasksum = 
          modMask(mask[i-1][j+1]) + modMask(mask[i][j+1]) + modMask(mask[i+1][j+1]) +      
          modMask(mask[i-1][j])   +                       + modMask(mask[i+1][j])   +
          modMask(mask[i-1][j-1]) + modMask(mask[i][j-1]) + modMask(mask[i+1][j-1]);
        if (neighmasksum <= (7*MASK_SHEET + MASK_DRAGGING)) { 
          mask[i][j] = MASK_SHEET;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  // communicate, because computing eff visc requires neighbor vals of mask
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  return 0;
}



PetscErrorCode IceDragYieldModel::writeMatlabDragFile(const char *basename) {
  PetscErrorCode  ierr;
  PetscViewer     viewer;
  char b[PETSC_MAX_PATH_LEN];
  char mf[PETSC_MAX_PATH_LEN];  // Matlab format

  strncpy(b, basename, PETSC_MAX_PATH_LEN-4);
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", b, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  strcpy(mf, b);
  strcat(mf, ".m");
  ierr = verbPrintf(2,grid.com, "writing variables beta, betax, betay, balvel to file %s ...", mf); CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(grid.com, mf, &viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\ndisp('IceDragYieldModel output:')\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nclear year x y xx yy beta betax betay dragxy balvel\n");
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"year=%10.6f;\n",grid.p->year);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "x=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",
                                grid.p->Lx,grid.p->dx,grid.p->Lx);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "y=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",
                                grid.p->Ly,grid.p->dy,grid.p->Ly);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n\n");  CHKERRQ(ierr);

  ierr=PetscObjectSetName((PetscObject) g2,"balvel"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vmagbalvel, g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nbalvel = reshape(balvel,%d,%d);\n\n",
                                grid.p->Mx,grid.p->My);  CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) dragxy,"dragxy"); CHKERRQ(ierr);
  ierr = VecView(dragxy, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"betax = dragxy(1:2:%d);\n",(2*grid.p->Mx*grid.p->My)-1);
        CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nbetax = reshape(betax,%d,%d);\n\n",
                                grid.p->Mx,grid.p->My);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"betay = dragxy(2:2:%d);\n",2*grid.p->Mx*grid.p->My);
        CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nbetay = reshape(betay,%d,%d);\n\n",
                                grid.p->Mx,grid.p->My);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer, "beta = (max(betax,0.0) + max(betay,0.0))/2;\n"); CHKERRQ(ierr);

  ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);
  return 0;
}


PetscErrorCode IceDragYieldModel::computeScalarBeta() {
  PetscErrorCode  ierr;
  PetscTruth      useMacayealVelocitySAVE;
  PetscScalar     muSlidingSAVE, enhancementFactorSAVE;

  useMacayealVelocitySAVE = useMacayealVelocity;
  muSlidingSAVE = muSliding;
  enhancementFactorSAVE = enhancementFactor;

    ierr = verbPrintf(2,grid.com, 
          "computing velocities (with MacAyeal; including eff viscosity iteration) ...\n");
          CHKERRQ(ierr);
    setUseMacayealVelocity(PETSC_TRUE);
    ierr = velocity(true); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
    ierr = updateViewers(); CHKERRQ(ierr);

    ierr = verbPrintf(2,grid.com, "computing resulting effective viscosity at all grounded points ...");
          CHKERRQ(ierr);
    ierr = getEffectiveViscosity(); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
    ierr = updateViewers(); CHKERRQ(ierr);

    ierr = verbPrintf(2,grid.com, "computing deformational velocities (w/o MacAyeal) ...");
          CHKERRQ(ierr);
    setUseMacayealVelocity(PETSC_FALSE);
    setMuSliding(0.0);  // for deformational, just assume frozen bed
    setEnhancementFactor(0.8);  //  reduce amount of deformation to ascribe more
                                //  of mass-balance velocities to sliding
    ierr = velocity(true); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
    ierr = updateViewers(); CHKERRQ(ierr);

    ierr = verbPrintf(2,grid.com, 
          "computing (scalar) beta by subtracting deformation from balance and using MacAyeal eqns ...");
          CHKERRQ(ierr);
    ierr = readBalVels(); CHKERRQ(ierr);
    ierr = computeDragFromBalanceVelocity(); CHKERRQ(ierr);
    ierr = moveDragxytoDAbetas(); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
    ierr = updateViewers(); CHKERRQ(ierr);

//    ierr = verbPrintf(2,grid.com, "saving Matlab file ..."); CHKERRQ(ierr);
//    ierr = writeMatlabDragFile("dragfile"); CHKERRQ(ierr);
//    ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, "saving beta in NetCDF file ..."); CHKERRQ(ierr);
    ierr = writeBetaNC(); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);

    ierr = verbPrintf(2,grid.com, "updating mask using beta ..."); CHKERRQ(ierr);
    ierr = updateMaskFromBeta(); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
    ierr = updateViewers(); CHKERRQ(ierr);
    
    setUseMacayealVelocity(useMacayealVelocitySAVE);
    setMuSliding(muSlidingSAVE);
    setEnhancementFactor(enhancementFactorSAVE);
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

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    IceGrid g(com, rank, size);
    PetscInt   flowlawNumber = 0;  // use Paterson-Budd by default
    IceType*   ice;

    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1,com, "PANT (PISM Antarctic IceDragYieldModel) run mode\n"); CHKERRQ(ierr);
    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);
    IceDragYieldModel m(g, *ice);
    ierr = m.setFromOptions(); CHKERRQ(ierr);
    ierr = m.initFromOptions(); CHKERRQ(ierr);
    ierr = m.setSoundingFromOptions(); CHKERRQ(ierr);

    ierr = m.dragYieldOptions(); CHKERRQ(ierr);
    ierr = m.createDragYieldVecs(); CHKERRQ(ierr);
    
    // will allow choice of scalarBeta or vectorBeta or Plastic
    ierr = m.computeScalarBeta(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "running (as PISM) ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);

    ierr = m.writeFiles("unnamed"); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, " ... done.\n"); CHKERRQ(ierr);
    
    ierr = verbPrintf(2,com, "destroying IceDragYieldModel Vecs ..."); CHKERRQ(ierr);
    ierr = m.destroyDragYieldVecs(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, " done\n"); CHKERRQ(ierr);
    
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

