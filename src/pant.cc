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

     pant -if ant40km.nc -mv -gk -e 1.2 -y 1 -verbose -betaxy -balvel init.nc -obasal betaxymap

[CONCRETELY:
     pant -if ant150k_40kmSMMV.nc -mv -gk -e 1.2 -ksp_rtol 1e-6 -ocean_kill -verbose \   
         -betaxy -balvel init.nc -obasal betaxymap -d Bncm -y 1
]

Here ant40km.nc contains a saved model state with a vaguely 
reasonable temperature field.  The balance velocities are read in from init.nc.
Then the deformation (SIA) velocity is computed and subtracted from
balance velocities to get a sliding velocity field everywhere ice is grounded.
The MacAyeal equations are (i.e. the MacAyeal nonlinear operators) applied to 
determine the vector (beta_x,beta_y).  These are written to a NetCDF file betaxymap.nc.
The mask is updated based on the value of scalar beta.  Then a 1 year run is done 
using the vector beta values; the model should be nearly in steady state 
(although it will drift away as the temperature changes).

2. For distributed *scalar* beta computed from balance velocities:

     pant -if ant40km.nc -mv -gk -e 1.2 -y 1 -verbose -beta -balvel init.nc -obasal betamap

[CONCRETELY:
     pant -if ant150k_40kmSMMV.nc -mv -gk -e 1.2 -ksp_rtol 1e-6 -ocean_kill -verbose \   
         -beta -balvel init.nc -obasal betamap -d Bncm -y 1
]

As in 1. except that once the vector (beta_x,beta_y) is generated we combine: beta
= (u^2 beta_x + v^2 beta_y)/(u^2+v^2) with max of BETAMAX = 1.0e13.  beta is written to 
a NetCDF file betamap.nc.  Then a 1 year run is done using the distributed beta; 
the model should less in steady state than in 1.

3. For Schoof-type solution of plastic till problem, with yield stress computed from 
thermal model and estimate of basal till pore water pressure:

     pant -if ant40km.nc -mv -gk -e 1.2 -y 1 -verbose -super -plastic -obasal taucmap

Here the plastic till Schoof version of the MacAyeal equations are solved from the saved model
state ant40km.nc.  In particular, the vHmelt values and the basal temps (vT) from
ant40km.nc are used to compute vtauc.  This map is saved in taucmap.nc.  Then 1 years
are run.  During this run the velocities computed from SIA are added to those from MacAyeal-Schoof
in the grounded ice.  The ice shelves have their usual scheme.
*/

#include <cstring>
#include <petscbag.h>
#include <netcdf.h>

#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"
#include "nc_util.hh"

class IceDragYieldModel : public IceModel {
public:
  PetscTruth     doPlastic;
  IceDragYieldModel(IceGrid &g, IceType &i);
  PetscErrorCode dragYieldInitFromOptions();
  PetscErrorCode dragYieldFinalize();
  PetscErrorCode writeBetaTaucNCFile();
  PetscErrorCode setupPlasticTauc();  // call to start plastic case (doPlastic == PETSC_TRUE)
  PetscErrorCode computeBeta();       // call to start linear case (doPlastic == PETSC_FALSE)
 
private:
  // note vWork2d[0:3] get used by methods of IceDragYieldModel
  Vec            vbetaxORbeta, vbetayORtauc;  // copies; see overloading comment above
  Vec            vNu[2];
  Vec 	         dragxy;    // this is not da2 sized; it's a (2 Mx My) x 1 column
  PetscTruth     doBetaxy;  // when doBetaxy = TRUE there is overloading of Vecs,
                            // namely vbeta stores beta_x while vtauc stores beta_y;
                            // this affects meaning of viewers "-d B" and "-d C"
  
  // for both plastic and linear (distributed beta) cases
  char           betataucFileName[PETSC_MAX_PATH_LEN];
  PetscErrorCode createDragYieldVecs();
  PetscErrorCode destroyDragYieldVecs();
  virtual PetscScalar basalDragx(PetscScalar**, PetscScalar**,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  virtual PetscScalar basalDragy(PetscScalar**, PetscScalar**,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  
  // for plastic case we insert basal water dependent computation of vtauc at start of each time step
  virtual PetscErrorCode additionalAtStartTimestep();
  
  // for distributed scalar beta or vector (beta_x,beta_y) basal (force linear in velocity)
  PetscTruth     balvelFileSet;
  char           balvelFileName[PETSC_MAX_PATH_LEN];
  Vec            vmagbalvel;
  PetscErrorCode readBalVels();
  PetscErrorCode getEffectiveViscosityAllGrounded();
  PetscErrorCode computeDragFromBalanceVelocity();
  PetscErrorCode updateMaskFromBeta();
  PetscErrorCode moveDragxytoDAbetas();
  PetscErrorCode replaceBadwMissing(Vec,bool);
};


// these are limits on computed and reported scalar, and components of vector, beta
// note units of BETAMAX/MIN are Pa s m^-1; BETAMAX is max allowed drag coeff beta
// BETAMAX is also missing_value in NetCDF output
#define BETAMAX 1.0e13
#define BETAMIN 1.0e5
// if mask[i][j] in these limits, and grounded, then marked as MASK_DRAGGING
#define MASK_BETA_RANGE_LOW   (5.0 * BETAMIN)
#define MASK_BETA_RANGE_HIGH  (0.01 * BETAMAX)
#define MASK_BETA_MAXRATIO    4.0

IceDragYieldModel::IceDragYieldModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {
  // do nothing (except pass g,i to constructor for IceModel)
}


PetscErrorCode IceDragYieldModel::createDragYieldVecs() {
  PetscErrorCode ierr;
  
  ierr = VecDuplicate(vh, &vbetaxORbeta); CHKERRQ(ierr);
  ierr = VecSet(vbetaxORbeta,BETAMAX); CHKERRQ(ierr);  // set to missing_value for now
  ierr = VecDuplicate(vh, &vbetayORtauc); CHKERRQ(ierr);
  ierr = VecSet(vbetayORtauc,BETAMAX); CHKERRQ(ierr);
  const PetscInt M = 2 * grid.p->Mx * grid.p->My;
  ierr = VecCreateMPI(grid.com, PETSC_DECIDE, M, &dragxy); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragYieldModel::destroyDragYieldVecs() {
  PetscErrorCode  ierr;

  ierr = VecDestroy(vbetaxORbeta); CHKERRQ(ierr);
  ierr = VecDestroy(vbetayORtauc); CHKERRQ(ierr);
  ierr = VecDestroy(dragxy); CHKERRQ(ierr);
  // note vmagbalvel, vNu[0], vNu[1] are just IceModel work Vecs
  return 0;
}


PetscErrorCode IceDragYieldModel::dragYieldInitFromOptions() {
  PetscErrorCode  ierr;
  PetscTruth      obasalSet, betaSet;

  ierr = PetscOptionsHasName(PETSC_NULL, "-beta", &betaSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-betaxy", &doBetaxy); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-plastic", &doPlastic); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-balvel", balvelFileName,
                               PETSC_MAX_PATH_LEN, &balvelFileSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-obasal", betataucFileName,
                               PETSC_MAX_PATH_LEN, &obasalSet); CHKERRQ(ierr);
  strcat(betataucFileName, ".nc");

  // throw errors if options conflict or don't include filenames
  if (doPlastic == PETSC_TRUE) {
    if (betaSet == PETSC_TRUE) {
      SETERRQ(1,"option conflict in pant; both -plastic and -beta set");
    }
    if (doBetaxy == PETSC_TRUE) {
      SETERRQ(2,"option conflict in pant; both -plastic and -betaxy set");
    }
    if (balvelFileSet == PETSC_TRUE) {
      verbPrintf(1,grid.com,
         "warning: both -plastic and -balvel set; ignoring balance velocities file");
    }
  } else {
    if (balvelFileSet == PETSC_FALSE) {
      SETERRQ(4,"pant requires source for balance velocities if -plastic is not set; use -balvel foo.nc");
    }
    if (betaSet == PETSC_TRUE) {
      if (doBetaxy == PETSC_TRUE) {
        SETERRQ(5,"option conflict in pant; both -beta and -betaxy set");
      }
    }
  }

  ierr = createDragYieldVecs(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragYieldModel::dragYieldFinalize() {
  PetscErrorCode  ierr;

  if (doPlastic == PETSC_TRUE) { // in this case we write the final tauc
    ierr = verbPrintf(2,grid.com, "saving tauc, beta, betax, betay in NetCDF file ..."); CHKERRQ(ierr);
    ierr = writeBetaTaucNCFile(); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
  }
  ierr = destroyDragYieldVecs(); CHKERRQ(ierr);
  return 0;
}


PetscScalar IceDragYieldModel::basalDragx(PetscScalar **betaORbetax, PetscScalar **taucORbetay,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  if (doBetaxy == PETSC_TRUE) {
    return betaORbetax[i][j]; // it is betax here
  } else { // In both plastic and scalar beta case, just use drag method of BasalType
    // class (see materials.hh).  Note this just returns beta[i][j] if doPlastic == FALSE 
    // while it returns tauc[i][j] over (regularized) speed if doPlastic == TRUE.
    return basal->drag(betaORbetax[i][j], taucORbetay[i][j], u[i][j], v[i][j]);  // beta and tauc here
  }
}


PetscScalar IceDragYieldModel::basalDragy(PetscScalar **betaORbetax, PetscScalar **taucORbetay,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  if (doBetaxy == PETSC_TRUE) {
    return taucORbetay[i][j]; // it is betay here
  } else { // SEE COMMENT IN basalDragx
    return basal->drag(betaORbetax[i][j], taucORbetay[i][j], u[i][j], v[i][j]);  // beta and tauc here
  }
}


PetscErrorCode nc_check(int stat) {
  if (stat)
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
  return 0;
}


PetscErrorCode IceDragYieldModel::writeBetaTaucNCFile() {
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
   int tauc_id;

#  define RANK_x 1
#  define RANK_y 1
#  define RANK_t 1
#  define RANK_beta 3
#  define RANK_betax 3
#  define RANK_betay 3
#  define RANK_tauc 3

   /* variable shapes */
   int x_dims[RANK_x];
   int y_dims[RANK_y];
   int t_dims[RANK_t];
   int beta_dims[RANK_beta];
   int betax_dims[RANK_betax];
   int betay_dims[RANK_betay];
   int tauc_dims[RANK_tauc];

   /* enter define mode */
if (grid.rank == 0) {
   stat = nc_create(betataucFileName, NC_CLOBBER|NC_64BIT_OFFSET, &ncid);
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

   tauc_dims[0] = t_dim;
   tauc_dims[1] = x_dim;
   tauc_dims[2] = y_dim;
   stat = nc_def_var(ncid, "tauc", NC_FLOAT, RANK_tauc, tauc_dims, &tauc_id);
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

   float miss = BETAMAX;

   stat = nc_put_att_text(ncid, beta_id, "long_name", 30, "ice_basal_friction_coefficient");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, beta_id, "units", 8, "Pa s m-1");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_float(ncid, beta_id, "missing_value", NC_FLOAT, 1, &miss);
   stat = nc_put_att_text(ncid, betax_id, "long_name", 42, "ice_basal_friction_coefficient_x_component");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, betax_id, "units", 8, "Pa s m-1");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_float(ncid, betax_id, "missing_value", NC_FLOAT, 1, &miss);
   stat = nc_put_att_text(ncid, betay_id, "long_name", 42, "ice_basal_friction_coefficient_y_component");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, betay_id, "units", 8, "Pa s m-1");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_float(ncid, betay_id, "missing_value", NC_FLOAT, 1, &miss);
   stat = nc_put_att_text(ncid, tauc_id, "long_name", 27, "ice_basal_till_yield_stress");
   CHKERRQ(nc_check(stat));
   stat = nc_put_att_text(ncid, tauc_id, "units", 2, "Pa");
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
  Vec  to_put = vWork2d[3];
  if (doPlastic == PETSC_TRUE) {  
    ierr = VecSet(to_put,BETAMAX); CHKERRQ(ierr);
  } else {
    ierr = VecCopy((doBetaxy == PETSC_TRUE) ? vbetaxORbeta : vbeta,to_put); CHKERRQ(ierr);
    ierr = replaceBadwMissing(to_put,true); CHKERRQ(ierr);
  }
  ierr = put_local_var(&grid, ncid, beta_id, NC_FLOAT, grid.da2, to_put,
  		       g2, s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  if (doPlastic == PETSC_TRUE) {
    ierr = VecSet(to_put,BETAMAX); CHKERRQ(ierr);
  } else {
    ierr = VecCopy((doBetaxy == PETSC_TRUE) ? vbeta : vbetaxORbeta,to_put); CHKERRQ(ierr);
    ierr = replaceBadwMissing(to_put,true); CHKERRQ(ierr);
  }
  ierr = put_local_var(&grid, ncid, betax_id, NC_FLOAT, grid.da2, to_put,
  		       g2, s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  if (doPlastic == PETSC_TRUE) {
    ierr = VecSet(to_put,BETAMAX); CHKERRQ(ierr);
  } else {
    ierr = VecCopy((doBetaxy == PETSC_TRUE) ? vtauc : vbetayORtauc,to_put); CHKERRQ(ierr);
    ierr = replaceBadwMissing(to_put,true); CHKERRQ(ierr);
  }
  ierr = put_local_var(&grid, ncid, betay_id, NC_FLOAT, grid.da2, to_put,
  		       g2, s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  if (doPlastic == PETSC_FALSE) {
    ierr = VecSet(to_put,BETAMAX); CHKERRQ(ierr);
  } else {
    ierr = VecCopy(vtauc,to_put); CHKERRQ(ierr);
    ierr = replaceBadwMissing(to_put,false); CHKERRQ(ierr);
  }
  ierr = put_local_var(&grid, ncid, tauc_id, NC_FLOAT, grid.da2, to_put,
  		       g2, s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // We are done with these buffers
  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_close(ncid);
    CHKERRQ(nc_check(stat));
  }

  return 0;
}


PetscErrorCode IceDragYieldModel::replaceBadwMissing(Vec vmine, bool betaflag) {
  PetscErrorCode  ierr;
  PetscScalar     **vvv;

    ierr = DAVecGetArray(grid.da2, vmine, &vvv); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (vvv[i][j] <= -BETAMAX + 10.0)  vvv[i][j] = BETAMAX;
        if ((betaflag) && (vvv[i][j] <= BETAMIN + 1e-5))   vvv[i][j] = BETAMAX;
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vmine, &vvv); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragYieldModel::additionalAtStartTimestep() {
  PetscErrorCode  ierr;

  if (doPlastic == PETSC_TRUE) {

    // we implement formula (2.4) in C. Schoof 2006 "A variational approach to ice
    // stream flow", J. Fluid Mech. vol 556 pp 227--251:
    //  (2.4)   \tau_c = \mu (\rho g H - p_w)
    // we modify it by:
    //   1. adding a small till cohesion (see Paterson 3rd ed table 8.1)
    //   2. replacing   p_w --> \lambda p_w   where \lambda = Hmelt / DEFAULT_MAX_HMELT;
    //      thus 0 <= \lambda <= 1 and \lambda = 0 when bed is frozen 
    //   3. computing a porewater pressure p_w which is the max of (0.95 of overburden)
    //      and the porewater pressure computed by formula (4) in 
    //      C. Ritz et al 2001 J. G. R. vol 106 no D23 pp 31943--31964;
    //      the modification of this porewater pressure as in Lingle&Brown 1987 is not 
    //      implementable because the "elevation of the bed at the grounding line"
    //      is at an unknowable location (we are not doing a flow line model!)

    const PetscScalar plastic_till_c_0 = 20e3;  // Pa; 20kPa = 0.2 bar; cohesion of till
    const PetscScalar plastic_till_mu = 0.466307658156;  // = tan(25^o); till friction angle
    const PetscScalar porewater_gamma = 0.95; // max allowed fraction of overburden
    
    PetscScalar **mask, **tauc, **H, **Hmelt, **bed; 
    ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          tauc[i][j] = 0.0;  
        } else { // grounded
          mask[i][j] = MASK_DRAGGING;  // in Schoof model, everything is dragging, so force this
          const PetscScalar overburdenP = ice.rho * grav * H[i][j];
          const PetscScalar drivingP = - ocean.rho * grav * bed[i][j];
          const PetscScalar pw = PetscMax(porewater_gamma * overburdenP, drivingP);
          const PetscScalar lambda = Hmelt[i][j] / DEFAULT_MAX_HMELT;  // note Hmelt[i][j]=0 if frozen
          tauc[i][j] = plastic_till_c_0 + plastic_till_mu * (overburdenP - lambda * pw);
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
    // communication of mask deliberately skipped; will happen next time step if not sooner
  }
  // if (doPlastic == PETSC_FALSE) then do nothing extra for each timestep
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


PetscErrorCode IceDragYieldModel::getEffectiveViscosityAllGrounded() {
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

  
// compute drag from velocities and balance velocities
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

  // build approximation of MacAyeal-Morland equations (A x = rhs) at all grounded points
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
      // remove deformational from balance velocities; put in x
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
  ierr = DAVecGetArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vbeta : vbetaxORbeta, &betax); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vtauc : vbetayORtauc, &betay); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vbetaxORbeta : vbeta, &beta); CHKERRQ(ierr);
  const PetscInt     M = 2 * grid.p->My;
  const PetscScalar  epsuv = PetscSqr(1.0 / secpera);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
    
      // do bounds checking on betax, betay so that
      //     BETAMIN <= betax,betay <= BETAMAX; note BETAMAX is missing_value
      betax[i][j] = betaxylong[i*M + 2*j];
      // note that on NetCDF write, BETAMAX is made to be missing value
      if (isnan(betax[i][j]))  betax[i][j] = BETAMAX;
      if (isinf(betax[i][j]))  betax[i][j] = BETAMAX;
      if (betax[i][j] > BETAMAX)  betax[i][j] = BETAMAX;
      if (betax[i][j] < BETAMIN)  betax[i][j] = BETAMIN;
      betay[i][j] = betaxylong[i*M + 2*j+1];
      if (isnan(betay[i][j]))  betay[i][j] = BETAMAX;
      if (isinf(betay[i][j]))  betay[i][j] = BETAMAX;
      if (betay[i][j] > BETAMAX)  betay[i][j] = BETAMAX;
      if (betay[i][j] < BETAMIN)  betay[i][j] = BETAMIN;

      const PetscScalar usqr = PetscSqr(u[i][j])+epsuv, vsqr = PetscSqr(v[i][j])+epsuv; 
      beta[i][j] = (usqr * betax[i][j] + vsqr * betay[i][j]) / (usqr + vsqr);
      if (beta[i][j] > BETAMAX)  beta[i][j] = BETAMAX;  // presumably redundant; accounts for possible rounding
      if (beta[i][j] < BETAMIN)  beta[i][j] = BETAMIN;  
    }
  }
  ierr = VecRestoreArray(dragxyLoc, &betaxylong); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vbeta : vbetaxORbeta, &betax); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vtauc : vbetayORtauc, &betay); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vbetaxORbeta : vbeta, &beta); CHKERRQ(ierr);

  ierr = VecDestroy(dragxyLoc); CHKERRQ(ierr);
  ierr = VecScatterDestroy(dragxyScatterGlobalToLocal); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragYieldModel::updateMaskFromBeta() {
  PetscErrorCode  ierr;
  PetscScalar     **mask, **beta, **betax, **betay, **newmask;
  Vec             vnewMask;

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vbetaxORbeta : vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vbeta : vbetaxORbeta, &betax); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vtauc : vbetayORtauc, &betay); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) != MASK_FLOATING) {
        // only mark as an ice stream if computed beta was reasonable
        const PetscScalar ratio = betax[i][j]/betay[i][j];
        if ( (beta[i][j] >= MASK_BETA_RANGE_LOW) && (beta[i][j] <= MASK_BETA_RANGE_HIGH) 
             && (ratio <= MASK_BETA_MAXRATIO) && ((1.0/ratio) <= MASK_BETA_MAXRATIO)     ) {
          mask[i][j] = MASK_DRAGGING;
        } else {
          mask[i][j] = MASK_SHEET;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vbetaxORbeta : vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vbeta : vbetaxORbeta, &betax); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, (doBetaxy == PETSC_TRUE) ? vtauc : vbetayORtauc, &betay); CHKERRQ(ierr);

  // 1.  remove singleton and nearly singleton DRAGGING points, i.e. if 
  //     all box stencil neighbors are SHEET or if at most one is DRAGGING or FLOATING
  // 2.  make SHEET points with >=3 box stencil FLOATING neighbors *and* 
  //     >=2 box stencil dragging neighbors into DRAGGING (i.e. don't block streams as they hit floating)
  ierr = VecDuplicate(vMask,&vnewMask); CHKERRQ(ierr);
  ierr = VecCopy(vMask,vnewMask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vnewMask, &newmask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const int neighmasksum = 
        modMask(mask[i-1][j+1]) + modMask(mask[i][j+1]) + modMask(mask[i+1][j+1]) +      
        modMask(mask[i-1][j])   +                       + modMask(mask[i+1][j])   +
        modMask(mask[i-1][j-1]) + modMask(mask[i][j-1]) + modMask(mask[i+1][j-1]);
      if ((intMask(mask[i][j]) == MASK_DRAGGING) && (neighmasksum <= (7*MASK_SHEET + MASK_FLOATING))) { 
        newmask[i][j] = MASK_SHEET;
      } else if ( (intMask(mask[i][j]) == MASK_SHEET)
                  && (neighmasksum >= (3*MASK_FLOATING + 2*MASK_DRAGGING + 3*MASK_SHEET)) ) { 
        newmask[i][j] = MASK_DRAGGING;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vnewMask, &newmask); CHKERRQ(ierr);
  ierr = VecCopy(vnewMask,vMask); CHKERRQ(ierr);

  // 3.  one more round of strict singleton removal for DRAGGING points surrounded by SHEET
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vnewMask, &newmask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const int neighmasksum = 
        modMask(mask[i-1][j+1]) + modMask(mask[i][j+1]) + modMask(mask[i+1][j+1]) +      
        modMask(mask[i-1][j])   +                       + modMask(mask[i+1][j])   +
        modMask(mask[i-1][j-1]) + modMask(mask[i][j-1]) + modMask(mask[i+1][j-1]);
      if ((intMask(mask[i][j]) == MASK_DRAGGING) && (neighmasksum <= (8*MASK_SHEET))) { 
        newmask[i][j] = MASK_SHEET;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vnewMask, &newmask); CHKERRQ(ierr);
  ierr = VecCopy(vnewMask,vMask); CHKERRQ(ierr);
  ierr = VecDestroy(vnewMask); CHKERRQ(ierr);

  // communicate, because computing eff visc requires neighbor vals of mask
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragYieldModel::setupPlasticTauc() {
  // PetscErrorCode  ierr;

  delete basal;
  basal = new PlasticBasalType;
  
  // see IceDragYieldModel::additionalAtStartTimeStep() for computation of tau_c
  // from porewater pressure estimate (i.e. using vHmelt)

  SETERRQ(1,"PLASTIC CASE NOT IMPLEMENTED; quitting ...\n");
  return 0;
}


PetscErrorCode IceDragYieldModel::computeBeta() {
  PetscErrorCode  ierr;
  PetscTruth      useMacayealVelocitySAVE, doBetaxySAVE;
  PetscScalar     muSlidingSAVE;

  // save aspects of IceModel state before computing beta
  useMacayealVelocitySAVE = useMacayealVelocity;
  doBetaxySAVE = doBetaxy;
  muSlidingSAVE = muSliding;

  ierr = verbPrintf(2,grid.com, 
        "computing velocities (with MacAyeal; including eff viscosity iteration) ...\n");
        CHKERRQ(ierr);
  setUseMacayealVelocity(PETSC_TRUE);
  doBetaxy = PETSC_FALSE; // so  IceDragYieldModel::basalDrag[x|y]()  have usual meaning
                          // {i.e. same meaning as  IceModel::basalDrag[x|y]() }
  ierr = velocity(true); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "computing resulting effective viscosity at all grounded points ...");
        CHKERRQ(ierr);
  ierr = getEffectiveViscosityAllGrounded(); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "computing deformational velocities (w/o MacAyeal) ...");
        CHKERRQ(ierr);
  setUseMacayealVelocity(PETSC_FALSE);
  setMuSliding(0.0);  // for deformational, just assume frozen bed
//  setEnhancementFactor(0.8);  //  reduce amount of deformation to ascribe more
//                              //  of mass-balance velocities to sliding
  ierr = velocity(true); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, 
        "computing (scalar and vector) beta by subtracting deformation from balance and using MacAyeal eqns ...");
        CHKERRQ(ierr);
  ierr = readBalVels(); CHKERRQ(ierr);
  ierr = computeDragFromBalanceVelocity(); CHKERRQ(ierr);

  // restore state before putting draxy into da2 Vecs (and writing nc file and updating mask ...) 
  setUseMacayealVelocity(useMacayealVelocitySAVE);
  setMuSliding(muSlidingSAVE);
  doBetaxy = doBetaxySAVE;

  ierr = moveDragxytoDAbetas(); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "saving beta, betax, betay, tauc in NetCDF file ..."); CHKERRQ(ierr);
  ierr = writeBetaTaucNCFile(); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "updating mask using beta ..."); CHKERRQ(ierr);
  ierr = updateMaskFromBeta(); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, " done \n"); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);

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

    // special procedures for derived class IceDragYieldModel:
    ierr = m.dragYieldInitFromOptions(); CHKERRQ(ierr);    
    if (m.doPlastic == PETSC_TRUE) {
      ierr = m.setupPlasticTauc(); CHKERRQ(ierr);
    } else {
      ierr = m.computeBeta(); CHKERRQ(ierr); // includes NetCDF write of basal fields
    }

    ierr = verbPrintf(2,com, "running (as PISM) ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);
    ierr = m.writeFiles("pant_unnamed"); CHKERRQ(ierr);
    ierr = m.dragYieldFinalize(); CHKERRQ(ierr); // includes NetCDF write of basal fields in plastic case
    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

