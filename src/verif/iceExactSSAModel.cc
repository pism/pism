// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "exactTestsIJ.h"
#include "exactTestM.h"
#include "iceExactSSAModel.hh"

const PetscScalar IceExactSSAModel::m_schoof = 10; // (pure number)
const PetscScalar IceExactSSAModel::L_schoof = 40e3; // meters
const PetscScalar IceExactSSAModel::aspect_schoof = 0.05; // (pure)
const PetscScalar IceExactSSAModel::H0_schoof = aspect_schoof * L_schoof; 
                                       // = 2000 m THICKNESS
const PetscScalar IceExactSSAModel::B_schoof = 3.7e8; // Pa s^{1/3}; hardness 
                                       // given on p. 239 of Schoof; why so big?
const PetscScalar IceExactSSAModel::p_schoof = 4.0/3.0; // = 1 + 1/n

const PetscScalar IceExactSSAModel::LforJ = 300.0e3; // 300 km half-width
const PetscScalar IceExactSSAModel::LforM = 750.0e3; // 750 km half-width


IceExactSSAModel::IceExactSSAModel(IceGrid &g, char mytest) : IceModel(g) {
  test = mytest;
  
  config.set("max_iterations_ssa", 500);
  config.set_flag("use_ssa_velocity",         true);
  config.set_flag("use_constant_nuh_for_ssa", false);
  config.set_flag("do_plastic_till",          true); // correct for I, irrelevant for J and M
  config.set_flag("do_superpose",             false);
}

PetscErrorCode IceExactSSAModel::init_physics() {
  PetscErrorCode ierr;
  if (test == 'I')
    iceFactory.setType(ICE_CUSTOM);
  else
    iceFactory.setType(ICE_ARR);

  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  // If the user left things alone, we'll have a CustomGlenIce
  CustomGlenIce *cust = dynamic_cast<CustomGlenIce*>(ice);
  if (!cust) {
    ierr = verbPrintf(2,grid.com,"Warning, custom ice not in use, reported errors will not be correct\n",test); CHKERRQ(ierr);
  } else {
    // Use Schoof's parameter
    ierr = cust->setHardness(B_schoof);CHKERRQ(ierr);
  }

  // If the user changes settings with -ice_custom_XXX, they asked for it.  If you don't want to allow this, disable the
  // line below.
  ierr = ice->setFromOptions();CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactSSAModel::setFromOptions() {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,grid.com,"initializing Test %c ... \n",test); CHKERRQ(ierr);

  // input file not allowed
  ierr = stop_if_set(grid.com, "-i"); CHKERRQ(ierr);
  ierr = stop_if_set(grid.com, "-boot_from"); CHKERRQ(ierr);

  ierr = IceModel::setFromOptions();CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceExactSSAModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  if (test == 'J') {
    ierr = vNuForJ[0].create(grid, "vNuForJ[0]", true); CHKERRQ(ierr);
    ierr = vNuForJ[1].create(grid, "vNuForJ[1]", true); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceExactSSAModel::misc_setup() {
  PetscErrorCode ierr;

  ierr = IceModel::misc_setup();


  switch (test) {
  case 'I':
      // fill vtauc with values for Schoof solution:
      ierr = taucSetI(); CHKERRQ(ierr);
      // set up remaining stuff:
      ierr = setInitStateAndBoundaryVelsI(); CHKERRQ(ierr);
      // set flags, parameters affecting solve of stream equations
      // so periodic grid works although h(-Lx,y) != h(Lx,y):
      config.set_flag("compute_surf_grad_inward_ssa", true);
      config.set("epsilon_ssa", 0.0);  // don't use this lower bound
    break;
  case 'J':
    ierr = setInitStateJ(); CHKERRQ(ierr);
    config.set_flag("is_dry_simulation", false);
    leaveNuHAloneSSA = PETSC_TRUE; // will use already-computed nu instead of updating
    config.set_flag("compute_surf_grad_inward_ssa", false);
    config.set("epsilon_ssa", 0.0);  // don't use this lower bound
    break;
  case 'M':
    ierr = setInitStateM(); CHKERRQ(ierr);
    config.set_flag("is_dry_simulation", false);
    config.set_flag("compute_surf_grad_inward_ssa", false);

    ierr = ice->printInfo(3);CHKERRQ(ierr);
    // EXPERIMENT WITH STRENGTH BEYOND CALVING FRONT:
    ierr = ssaStrengthExtend.set_notional_strength(6.5e+16);CHKERRQ(ierr); // about optimal; compare 4.74340e+15 usual
    ierr = verbPrintf(3,grid.com,
		      "IceExactSSAModel::misc_setup, for test M:\n"
		      "  use_constant_nuh_for_ssa=%d\n",
		      config.get_flag("use_constant_nuh_for_ssa")); CHKERRQ(ierr);
  }

  // Communicate so that we can differentiate surface, and to set boundary conditions
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[0].endGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[1].endGhostComm(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceExactSSAModel::set_grid_defaults() {
  PetscErrorCode ierr;
  const PetscScalar testI_Ly = 120.0e3,
    testI_Lx = PetscMax(60.0e3, ((grid.Mx - 1) / 2) * (2.0 * testI_Ly / (grid.My - 1)) );

  switch (test) {
  case 'I':
    // set up X km by 240 km by 3000 m grid; note 240km = 6L where L = L_schoof
    // find dimension X so that pixels are square but for some unknown reason 
    //    there are convergence problems when x dimension Lx is significantly
    //    smaller than y dimension Ly
    // note Mx,My may have been set by options, but not others
    grid.Lx = testI_Lx;
    grid.Ly = testI_Ly;
    grid.Lz = 3000.0;
    break;
  case 'J':
    grid.Lx = grid.Ly = LforJ;
    // grid is truly periodic both in x and in y directions
    grid.periodicity = XY_PERIODIC;
    break;
  case 'M':
    grid.Lx = grid.Ly = LforM;
    grid.Lz = 1000.0;
    break;
  default:
    ierr = PetscPrintf(grid.com, "Only tests I, J, M currently supported in IceExactSSAModel.");
    CHKERRQ(ierr);
    PetscEnd();
  }

  return 0;
}

PetscErrorCode IceExactSSAModel::set_vars_from_options() {
  PetscErrorCode ierr;

  // We need a pointer to surface temp from PISMAtmosphereCoupler atmosPCC*
  IceModelVec2  *pccTs;
  ierr = atmosPCC->updateSurfTempAndProvide(grid.year, 0.0, // year and dt are irrelevant here 
					    pccTs); CHKERRQ(ierr);  

  // fill in temperature and age; not critical
  const PetscScalar T0 = 263.15;  // completely arbitrary
  ierr = pccTs->set(T0); CHKERRQ(ierr);
  ierr =  T3.set(T0); CHKERRQ(ierr);
  ierr = Tb3.set(T0); CHKERRQ(ierr);

  ierr = tau3.set(0.0); CHKERRQ(ierr);  // age (not yield stress)

  // set initial velocities (for start of iteration)
  ierr = vubar.set(0.0); CHKERRQ(ierr);
  ierr = vvbar.set(0.0); CHKERRQ(ierr);
  ierr = vuvbar[0].set(0.0); CHKERRQ(ierr);
  ierr = vuvbar[1].set(0.0); CHKERRQ(ierr);
  // clear 3D and basal velocities too
  ierr = u3.set(0.0); CHKERRQ(ierr);
  ierr = v3.set(0.0); CHKERRQ(ierr);
  ierr = w3.set(0.0); CHKERRQ(ierr);
  ierr = vub.set(0.0); CHKERRQ(ierr);
  ierr = vvb.set(0.0); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceExactSSAModel::taucSetI() {
  PetscErrorCode ierr;
  PetscScalar **tauc;

  ierr = vtauc.get_array(tauc); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt jfrom0 = j - (grid.My - 1)/2;
      const PetscScalar y = grid.dy * jfrom0;
      const PetscScalar theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
      const PetscScalar f = ice->rho * earth_grav * H0_schoof * tan(theta);
      tauc[i][j] = f * pow(PetscAbs(y / L_schoof), m_schoof);
    }
  }
  ierr = vtauc.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::setInitStateAndBoundaryVelsI() {
  PetscErrorCode ierr;
  PetscScalar    **uvbar[2], **mask, **h, **bed;
  
  ierr = vMask.set(MASK_DRAGGING); CHKERRQ(ierr);
  ierr = vH.set(H0_schoof); CHKERRQ(ierr);

  // set h, bed everywhere
  // on edges y = +- 3 L_schoof, set velocity and make mask=SHEET
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);    
  ierr = vbed.get_array(bed); CHKERRQ(ierr);
  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk, myu, myv;
      const PetscInt ifrom0 = i - (grid.Mx - 1)/2,
                     jfrom0 = j - (grid.My - 1)/2;
      const PetscScalar myx = grid.dx * ifrom0, myy = grid.dy * jfrom0;
      // eval exact solution; will only use exact vels if at edge
      exactI(m_schoof, myx, myy, &(bed[i][j]), &junk, &myu, &myv); 
      h[i][j] = bed[i][j] + H0_schoof;
      bool edge = ( (j == 0) || (j == grid.My - 1) );
      if (edge) {
        // set boundary condition which will apply to finite difference system:
        // staggered grid velocities at MASK_SHEET points at edges of grid
        mask[i][j] = MASK_SHEET;
        uvbar[0][i-1][j] = myu;
        uvbar[0][i][j] = myu;    // so average onto regular grid point (i,j) has u=myu
        uvbar[1][i][j-1] = myv;
        uvbar[1][i][j] = myv;    // so average onto regular grid point (i,j) has v=myv
      }
    }
  }  
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);    
  ierr = vbed.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::setInitStateJ() {
  PetscErrorCode ierr;
  PetscScalar    **H, **h, **mask, **uvbar[2];

  ierr = vbed.set(-5000.0); CHKERRQ(ierr); // assures shelf is floating
  ierr = vMask.set(MASK_FLOATING); CHKERRQ(ierr);

  /* use Ritz et al (2001) value of 30 MPa yr for typical vertically-averaged viscosity */
  const PetscScalar nu0 = 30.0 * 1.0e6 * secpera; /* = 9.45e14 Pa s */
  const PetscScalar H0 = 500.0;       /* 500 m typical thickness */
  ierr = vNuForJ[0].set(nu0 * H0); CHKERRQ(ierr);
  ierr = vNuForJ[1].set(nu0 * H0); CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);    
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk1, myu, myv;
      const PetscInt ifrom0 = i - (grid.Mx)/2,
                     jfrom0 = j - (grid.My)/2;
      const PetscScalar myx = grid.dx * ifrom0, myy = grid.dy * jfrom0;
      // set H,h on regular grid
      ierr = exactJ(myx, myy, &H[i][j], &junk1, &myu, &myv); CHKERRQ(ierr);
      h[i][j] = (1.0 - ice->rho / ocean.rho) * H[i][j];
      // special case at center point: here we indirectly set ubar,vbar 
      // at (i,j) by marking this grid point as SHEET and setting staggered-grid
      // version of ubar approriately; the average done in assembleSSARhs()
      // for SHEET points puts our value of velocity onto regular grid point (i,j)
      // and makes ubar=myu, vbar=myv there
      if ( (i == (grid.Mx)/2) && (j == (grid.My)/2) ) {
        mask[i][j] = MASK_SHEET;
        uvbar[0][i-1][j] = myu;   uvbar[0][i][j] = myu;
        uvbar[1][i][j-1] = myv;   uvbar[1][i][j] = myv;
        
        ierr = u3.begin_access(); CHKERRQ(ierr);
        ierr = v3.begin_access(); CHKERRQ(ierr);
        ierr = u3.setColumn(i,j,myu); CHKERRQ(ierr);
        ierr = v3.setColumn(i,j,myv); CHKERRQ(ierr);
        ierr = u3.end_access(); CHKERRQ(ierr);
        ierr = v3.end_access(); CHKERRQ(ierr);
      }
    }
  }  
  ierr = vh.end_access(); CHKERRQ(ierr);    
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::setInitStateM() {
  PetscErrorCode ierr;
  PetscScalar    **H, **h, **bed, **mask, **ubar, **vbar,
                 **ub, **vb, **uvbar[2];

  const PetscScalar
            Rg = 300.0e3,     // radius for grounding line
            Rc = 600.0e3,     // radius for calving front
            H0 = 500.0,       // 500 m constant thickness
            hicepresent = (1.0 - (ice->rho / ocean.rho)) * H0,
            bedgrounded = - H0 * (ice->rho / ocean.rho);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);    
  ierr = vbed.get_array(bed); CHKERRQ(ierr);    
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vubar.get_array(ubar); CHKERRQ(ierr);
  ierr = vvbar.get_array(vbar); CHKERRQ(ierr);
  ierr = vub.get_array(ub); CHKERRQ(ierr);
  ierr = vvb.get_array(vb); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar alpha, Drr, myu, myv;
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      // evaluate exact solution (though velocity info discarded for shelf)
      ierr = exactM(r, &alpha, &Drr, 1.0e-12, 0.0, 1); // CHKERRQ(ierr);
      if (ierr != 0) {
        verbPrintf(1,grid.com,"exactM() not successful; gives alpha = %f (m/a), at\n",
                   alpha*secpera);
        verbPrintf(1,grid.com,"   i,j,xx,yy,r = %d,%d,%f,%f,%f\n", i,j,xx,yy,r);
      }
      if (r > 1.0) { // merely to avoid div by 0
        myu = alpha * (xx / r); myv = alpha * (yy / r);
      } else {
        myu = 0.0; myv = 0.0;
      }
      if (r < Rg) {
        // grounded case: set velocities (to give boundary condition for shelf)
        H[i][j] = H0;
        h[i][j] = hicepresent;
        bed[i][j] = bedgrounded;
        mask[i][j] = MASK_SHEET;
        // see IceModel::broadcastSSAVelocity() to see how velocity will be set
        //   for floating portion; do that for grounded now
        ubar[i][j] = myu; 
        vbar[i][j] = myv;
        ub[i][j] = myu;
        vb[i][j] = myv;
        ierr = u3.setColumn(i,j,myu); CHKERRQ(ierr);
        ierr = v3.setColumn(i,j,myv); CHKERRQ(ierr);
        // w3 is set by IceModel::vertVelocityFromIncompressibility()
        // see IceModel::assembleSSARhs() to see why pairs are set here:
        uvbar[0][i][j] = myu;  uvbar[0][i-1][j] = myu;
        uvbar[1][i][j] = myv;  uvbar[1][i][j-1] = myv;
      } else if (r <= Rc) {
        // ice shelf case
        H[i][j] = H0;
        h[i][j] = hicepresent;
        // drop away at 2 m per hor. km beyond grounding line:
        bed[i][j] = bedgrounded - (r - Rg) * 0.002;  
        mask[i][j] = MASK_FLOATING;
      } else { // r > Rc
        // ocean case
        H[i][j] = 0.0;
        h[i][j] = 0.0;
        // continue to drop away at 2 m per hor. km beyond grounding line:
        bed[i][j] = bedgrounded - (r - Rg) * 0.002;  
        mask[i][j] = MASK_FLOATING_OCEAN0;
      }
    }
  }  
  ierr = vh.end_access(); CHKERRQ(ierr);    
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);    
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vubar.end_access(); CHKERRQ(ierr);
  ierr = vvbar.end_access(); CHKERRQ(ierr);
  ierr = vub.end_access(); CHKERRQ(ierr);
  ierr = vvb.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::reportErrors() {
  PetscErrorCode  ierr;
  PetscScalar exactmaxu, maxvecerr = 0.0, avvecerr = 0.0, 
              avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  PetscScalar gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0,
              gmaxuerr = 0.0, gmaxverr = 0.0;
  PetscScalar **u, **v;

  if (config.get_flag("do_pseudo_plastic_till")) {
    ierr = verbPrintf(1,grid.com, 
          "WARNING: numerical errors not valid for pseudo-plastic till\n"); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com, 
          "NUMERICAL ERRORS in velocity relative to exact solution:\n"); CHKERRQ(ierr);

  ierr = vubar.get_array(u); CHKERRQ(ierr);
  ierr = vvbar.get_array(v); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2, uexact, vexact;
      PetscScalar myr,myx,myy;
      mapcoords(i,j,myx,myy,myr);
      // eval exact solution
      if (test == 'I') {
        exactI(m_schoof, myx, myy, &junk1, &junk2, &uexact, &vexact); 
      } else if (test == 'J') {
        const PetscInt ifrom0 = i - (grid.Mx)/2,
                       jfrom0 = j - (grid.My)/2;
        const PetscScalar myx = grid.dx * ifrom0, myy = grid.dy * jfrom0;
        exactJ(myx, myy, &junk1, &junk2, &uexact, &vexact);
      } else if (test == 'M') {
        PetscScalar alpha;
        ierr = exactM(myr, &alpha, &junk1, 1.0e-12, 0.0, 1); //CHKERRQ(ierr);
        if (ierr != 0) {
          verbPrintf(1,grid.com,"exactM() not successful; gives alpha = %f (m/a), at\n",
                     alpha*secpera);
          verbPrintf(1,grid.com,"   i,j,xx,yy,r = %d,%d,%f,%f,%f\n", i,j,myx,myy,myr);
        }
        if (myr > 1.0) {
          uexact = alpha * (myx / myr);
          vexact = alpha * (myy / myr);
        } else {
          uexact = 0.0; vexact = 0.0;
        }
      } else {
        SETERRQ(1,"only tests I,J,M have computable errors");
      }
      // compute maximum errors
      const PetscScalar uerr = PetscAbsReal(u[i][j] - uexact);
      const PetscScalar verr = PetscAbsReal(v[i][j] - vexact);
      avuerr = avuerr + uerr;      
      avverr = avverr + verr;      
      maxuerr = PetscMax(maxuerr,uerr);
      maxverr = PetscMax(maxverr,verr);
      const PetscScalar vecerr = sqrt(uerr * uerr + verr * verr);
      maxvecerr = PetscMax(maxvecerr,vecerr);
      avvecerr = avvecerr + vecerr;
    }
  }
  ierr = vubar.end_access(); CHKERRQ(ierr);
  ierr = vvbar.end_access(); CHKERRQ(ierr);
   
  ierr = PetscGlobalMax(&maxuerr, &gmaxuerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxverr, &gmaxverr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avuerr, &gavuerr, grid.com); CHKERRQ(ierr);
  gavuerr = gavuerr/(grid.Mx*grid.My);
  ierr = PetscGlobalSum(&avverr, &gavverr, grid.com); CHKERRQ(ierr);
  gavverr = gavverr/(grid.Mx*grid.My);
  ierr = PetscGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.Mx*grid.My);

  if (test == 'I') {
    PetscScalar junk1, junk2, junk3;
    exactI(m_schoof, 0.0, 0.0, &junk1, &junk2, &exactmaxu, &junk3);
    ierr = verbPrintf(1,grid.com, 
       "velocity  :  maxvector   prcntavvec      maxu       avu\n");
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
            "           %11.4f%13.5f%10.4f%10.4f\n", 
            gmaxvecerr*secpera, (gavvecerr/exactmaxu)*100.0,
            gmaxuerr*secpera, gavuerr*secpera); CHKERRQ(ierr);
  } else if ((test == 'J') || (test == 'M')) {
    if (test == 'J') {
      // following from "pismv -test J -Mx 601 -My 601 -Mz 3 -verbose -eo"
      exactmaxu = 181.366 / secpera;
    } else {
      // from simpleM with r=600.0:
      exactmaxu = 1180.79793 / secpera;
    }
    ierr = verbPrintf(1,grid.com, 
       "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
            "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n", 
            gmaxvecerr*secpera, (gavvecerr/exactmaxu)*100.0,
            gmaxuerr*secpera, gmaxverr*secpera, gavuerr*secpera, 
            gavverr*secpera); CHKERRQ(ierr);
  }


  if (test == 'I') {
    // also print max of u
    ierr = verbPrintf(3,grid.com, 
       "(exact maximum of u is %11.4f (m/a))\n",exactmaxu*secpera); CHKERRQ(ierr);
  }

  ierr = verbPrintf(1,grid.com, "NUM ERRORS DONE\n");  CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::fillFromExactSolution() {
  PetscErrorCode  ierr;
  PetscScalar **ubar, **vbar;

  ierr = vubar.get_array(ubar); CHKERRQ(ierr);
  ierr = vvbar.get_array(vbar); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2;
      PetscScalar myr,myx,myy;
      mapcoords(i,j,myx,myy,myr);
      if (test == 'I') {
        exactI(m_schoof, myx, myy, &junk1, &junk2, &ubar[i][j], &vbar[i][j]); 
      } else if (test == 'J') {
        const PetscInt ifrom0 = i - (grid.Mx)/2,
                       jfrom0 = j - (grid.My)/2;
        const PetscScalar myx = grid.dx * ifrom0, myy = grid.dy * jfrom0;
        ierr = exactJ(myx, myy, &junk1, &junk2, &ubar[i][j], &vbar[i][j]); CHKERRQ(ierr);
      } else if (test == 'M') {
        PetscScalar alpha;
        ierr = exactM(myr, &alpha, &junk1, 1.0e-12, 0.0, 1); //CHKERRQ(ierr);
        if (ierr != 0) {
          verbPrintf(1,grid.com,"exactM() not successful; gives alpha = %f (m/a), at\n",
                     alpha*secpera);
          verbPrintf(1,grid.com,"   i,j,xx,yy,r = %d,%d,%f,%f,%f\n", i,j,myx,myy,myr);
        }
        if (myr > 1.0) {
          ubar[i][j] = alpha * (myx / myr);
          vbar[i][j] = alpha * (myy / myr);
        } else {
          ubar[i][j] = 0.0;
          vbar[i][j] = 0.0;
        }
      } else {
        SETERRQ(1,"only tests I,J,M supported in IceExactSSAModel");
      }
      // now fill in velocity at depth
      for (PetscInt k=0; k<grid.Mz; k++) {
        ierr = u3.setColumn(i,j,ubar[i][j]); CHKERRQ(ierr);
        ierr = v3.setColumn(i,j,vbar[i][j]); CHKERRQ(ierr);
        // leave w alone for now; will be filled by call to broadcastSSAVelocity()
      }
    }
  }
  ierr = vubar.end_access(); CHKERRQ(ierr);
  ierr = vvbar.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::diagnosticRun() {
  PetscErrorCode  ierr;

  ierr = verbPrintf(2,grid.com, "running Test %c ...\n", test); CHKERRQ(ierr);

  // does user want to turn off actual numerical evolution (and save exact solution)?
  ierr = check_option("-eo", exactOnly); CHKERRQ(ierr);
  if (exactOnly == PETSC_TRUE) {
    ierr=verbPrintf(2,grid.com,"  EXACT SOLUTION ONLY, NO NUMERICAL SOLUTION\n"); CHKERRQ(ierr);
  }

  // print out some stats about input state
  ierr = summaryPrintLine(PETSC_TRUE,PETSC_TRUE, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
           CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep
  skipCountDown = 0;

  if (exactOnly == PETSC_TRUE) {
    // just fill with exact solution, including exact 3D hor velocities
    ierr = fillFromExactSolution(); CHKERRQ(ierr);
  } else { // numerically solve ice shelf/stream equations
    PetscInt numiter;
    ierr = initSSA(); CHKERRQ(ierr);
    if ((test == 'I') || (test == 'M')) {
      ierr = velocitySSA(&numiter); CHKERRQ(ierr);
    } else if (test == 'J') {
      // use locally allocated space for (computed) nu:
      PetscTruth     dosnes;
      ierr = check_option("-ssa_snes", dosnes); CHKERRQ(ierr);  
      if (dosnes == PETSC_TRUE) {
        ierr = velocitySSA_SNES(vNuForJ, &numiter); CHKERRQ(ierr); 
      } else {
        ierr = velocitySSA(vNuForJ, &numiter); CHKERRQ(ierr);
      }
      ierr = vNuForJ[0].destroy(); CHKERRQ(ierr); // immediately de-allocate
      ierr = vNuForJ[1].destroy(); CHKERRQ(ierr); // immediately de-allocate
    }
    // fill in 3D velocities (u,v,w)
    ierr = broadcastSSAVelocity(true); CHKERRQ(ierr);
  }

  // finally update w
  ierr = u3.beginGhostComm(); CHKERRQ(ierr);
  ierr = u3.endGhostComm(); CHKERRQ(ierr);
  ierr = v3.beginGhostComm(); CHKERRQ(ierr);
  ierr = v3.endGhostComm(); CHKERRQ(ierr);
  ierr = vertVelocityFromIncompressibility(); CHKERRQ(ierr);

  // report on result of computation (i.e. to standard out and to viewers)
  ierr = computeMax3DVelocities(); CHKERRQ(ierr); 
  ierr = summary(true,true); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);
  PetscInt    pause_time = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, PETSC_NULL); CHKERRQ(ierr);
  if (pause_time > 0) {
    ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }
  return 0;
}


void IceExactSSAModel::mapcoords(PetscInt i, PetscInt j,
                                 PetscScalar &x, PetscScalar &y, PetscScalar &r) {
  // compute x,y,r on grid from i,j
  PetscScalar ifrom0, jfrom0;

  ifrom0=static_cast<PetscScalar>(i)-static_cast<PetscScalar>(grid.Mx - 1)/2.0;
  jfrom0=static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
  x=grid.dx*ifrom0;
  y=grid.dy*jfrom0;
  r = sqrt(PetscSqr(x) + PetscSqr(y));
}

