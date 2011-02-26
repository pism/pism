// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "grid.hh"
#include "iceModel.hh"
#include "tests/exactTestsIJ.h"
#include "tests/exactTestM.h"
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


IceExactSSAModel::IceExactSSAModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides, char mytest)
  : IceModel(g, conf, conf_overrides) {
  test = mytest;
  
  exactOnly = PETSC_FALSE;
  dt = 0;

  config.set("max_iterations_ssa", 500);
  config.set_flag("use_ssa_velocity",         true);
  config.set_flag("use_constant_nuh_for_ssa", false);
  config.set_flag("use_ssa_when_grounded",    true); // correct for I, irrelevant for J and M
}

PetscErrorCode IceExactSSAModel::setFromOptions() {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,grid.com,"initializing Test %c ... \n",test); CHKERRQ(ierr);

  // input file not allowed
  ierr = stop_if_set(grid.com, "-i"); CHKERRQ(ierr);
  ierr = stop_if_set(grid.com, "-boot_file"); CHKERRQ(ierr);

  ierr = IceModel::setFromOptions();CHKERRQ(ierr);

  // can not be overridden 
  config.set_flag("do_cold_ice_methods", true);

  return 0;
}

PetscErrorCode IceExactSSAModel::init_physics() {
  PetscErrorCode ierr;

  iceFactory.setType(ICE_CUSTOM);

  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  // If the user left things alone, we'll have a CustomGlenIce
  CustomGlenIce *cust = dynamic_cast<CustomGlenIce*>(ice);
  if (!cust) {
    ierr = verbPrintf(2,grid.com,"Warning, custom ice not in use, reported errors will not be correct\n",test); CHKERRQ(ierr);
  } else {
    // Use Schoof's parameter
    cust->setHardness(B_schoof);
  }

  // If the user changes settings with -ice_custom_XXX, they asked for it.  If you don't want to allow this, disable the
  // line below.
  ierr = ice->setFromOptions();CHKERRQ(ierr);

  delete stress_balance;
  // the destructor os PISMStressBalance will take care of these:
  ssa = new SSAFD(grid, *basal, *ice, *EC, config);
  modifier = new SSBM_Trivial(grid, *ice, *EC, config);

  stress_balance = new PISMStressBalance(grid, ssa, modifier, config);

  ierr = stress_balance->init(variables); CHKERRQ(ierr);

  return 0;
}



PetscErrorCode IceExactSSAModel::misc_setup() {
  PetscErrorCode ierr;

  ierr = IceModel::misc_setup();

  switch (test) {
  case 'I':
    ierr = taucSetI(); CHKERRQ(ierr); // fill vtauc with values for Schoof soln
    ierr = setInitStateAndBoundaryVelsI(); CHKERRQ(ierr);
    // set flag so periodic grid works although h(-Lx,y) != h(Lx,y):
    config.set_flag("compute_surf_grad_inward_ssa", true);
    config.set("epsilon_ssa", 0.0);  // don't use this lower bound
    break;
  case 'J':
    ierr = setInitStateJ(); CHKERRQ(ierr);
    config.set_flag("is_dry_simulation", false);
    config.set_flag("compute_surf_grad_inward_ssa", false);
    config.set("epsilon_ssa", 0.0);  // don't use this lower bound
    // ensure that the strength extension's nu*H is used throughout:
    ierr = ssa->strength_extension->set_min_thickness(1000); CHKERRQ(ierr);
    break;
  case 'M':
    ierr = setInitStateM(); CHKERRQ(ierr);
    config.set_flag("is_dry_simulation", false);
    config.set_flag("compute_surf_grad_inward_ssa", false);

    // ierr = ice->printInfo(3);CHKERRQ(ierr);
    // EXPERIMENT WITH STRENGTH BEYOND CALVING FRONT:  correct value is unknown!!
    ierr = ssa->strength_extension->set_notional_strength(1.0e+15);CHKERRQ(ierr);
    ierr = verbPrintf(3,grid.com,
		      "IceExactSSAModel::misc_setup, for test M:\n"
		      "  use_constant_nuh_for_ssa=%d\n",
		      config.get_flag("use_constant_nuh_for_ssa")); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceExactSSAModel::set_grid_defaults() {
  PetscErrorCode ierr;
  const PetscScalar testI_Ly = 120.0e3,
    testI_Lx = PetscMax(60.0e3, ((grid.Mx - 1) / 2) * (2.0 * testI_Ly / (grid.My - 1)) );
  
  grid.start_year = grid.end_year = 0.0;

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
    PISMEnd();
  }

  return 0;
}

PetscErrorCode IceExactSSAModel::set_vars_from_options() {
  PetscErrorCode ierr;

  // fill in temperature and age; not critical
  const PetscScalar T0 = 263.15;  // completely arbitrary
  ierr = artm.set(T0); CHKERRQ(ierr);
  ierr =   T3.set(T0); CHKERRQ(ierr);
  ierr =  Tb3.set(T0); CHKERRQ(ierr);

  ierr = compute_enthalpy_cold(T3, Enth3); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceExactSSAModel::taucSetI() {
  PetscErrorCode ierr;
  PetscScalar **tauc;

  ierr = vtauc.get_array(tauc); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar y = grid.y[j];
      const PetscScalar theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
      const PetscScalar f = ice->rho * standard_gravity * H0_schoof * tan(theta);
      tauc[i][j] = f * pow(PetscAbs(y / L_schoof), m_schoof);
    }
  }
  ierr = vtauc.end_access(); CHKERRQ(ierr);
  // vtauc is global so no ghosts to update
  return 0;
}


PetscErrorCode IceExactSSAModel::setInitStateAndBoundaryVelsI() {
  PetscErrorCode ierr;
  PetscScalar    **h, **bed;
  IceModelVec2V &vel_bc = vWork2dV;
  
  ierr = vMask.set(MASK_DRAGGING_SHEET); CHKERRQ(ierr);
  ierr = vH.set(H0_schoof); CHKERRQ(ierr);

  // set h, bed everywhere
  // on edges y = +- 3 L_schoof, set velocity and make mask=SHEET
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);    
  ierr = vbed.get_array(bed); CHKERRQ(ierr);
  ierr = vel_bc.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk, myu, myv;
      const PetscScalar myx = grid.x[i], myy = grid.y[j];
      // eval exact solution; will only use exact vels if at edge
      exactI(m_schoof, myx, myy, &(bed[i][j]), &junk, &myu, &myv); 
      h[i][j] = bed[i][j] + H0_schoof;
      bool edge = ( (j == 0) || (j == grid.My - 1) );
      if (edge) {
        // set boundary condition which will apply to finite difference system:
        // staggered grid velocities at MASK_SHEET points at edges of grid
        vMask(i,j) = MASK_SHEET;
        vel_bc(i,j).u = myu;
        vel_bc(i,j).v = myv;
      }
    }
  }  
  ierr = vel_bc.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);    
  ierr = vbed.end_access(); CHKERRQ(ierr);

  ierr = stress_balance->set_boundary_conditions(vMask, vel_bc); CHKERRQ(ierr);

  // communicate what we have set
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);
  ierr = vbed.beginGhostComm(); CHKERRQ(ierr);
  ierr = vbed.endGhostComm(); CHKERRQ(ierr);
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::setInitStateJ() {
  PetscErrorCode ierr;
  IceModelVec2V &vel_bc = vWork2dV;

  ierr = vbed.set(-5000.0); CHKERRQ(ierr); // assures shelf is floating
  ierr = vMask.set(MASK_FLOATING); CHKERRQ(ierr);

  double ocean_rho = config.get("sea_water_density");

  /* use Ritz et al (2001) value of 30 MPa yr for typical vertically-averaged viscosity */
  const PetscScalar nu0 = 30.0 * 1.0e6 * secpera; /* = 9.45e14 Pa s */
  const PetscScalar H0 = 500.0;       /* 500 m typical thickness */

  ierr = ssa->strength_extension->set_notional_strength(nu0 * H0); CHKERRQ(ierr);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vh.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vel_bc.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk1, myu, myv;
      const PetscScalar myx = grid.x[i], myy = grid.y[j];
      // set H,h on regular grid
      ierr = exactJ(myx, myy, &vH(i,j), &junk1, &myu, &myv); CHKERRQ(ierr);
      vh(i,j) = (1.0 - ice->rho / ocean_rho) * vH(i,j);
      // special case at center point: here we set vel_bc at (i,j) by marking
      // this grid point as SHEET and setting ubar,vbar approriately
      if ( (i == (grid.Mx)/2) && (j == (grid.My)/2) ) {
        vMask(i,j) = MASK_SHEET;
        vel_bc(i,j).u = myu;
        vel_bc(i,j).v = myv;
      }
    }
  }  
  ierr = vh.end_access(); CHKERRQ(ierr);    
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vel_bc.end_access(); CHKERRQ(ierr);

  ierr = stress_balance->set_boundary_conditions(vMask, vel_bc); CHKERRQ(ierr);

  // communicate what we have set
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);
  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactSSAModel::setInitStateM() {
  PetscErrorCode ierr;
  PetscScalar    **H, **h, **bed, **mask;
  IceModelVec2V &vel_bc = vWork2dV;

  double ocean_rho = config.get("sea_water_density");

  const PetscScalar
            Rg = 300.0e3,     // radius for grounding line
            Rc = 600.0e3,     // radius for calving front
            H0 = 500.0,       // 500 m constant thickness
            hicepresent = (1.0 - (ice->rho / ocean_rho)) * H0,
            bedgrounded = - H0 * (ice->rho / ocean_rho);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);    
  ierr = vbed.get_array(bed); CHKERRQ(ierr);    
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vel_bc.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar alpha, Drr, myu, myv;
      PetscScalar xx = grid.x[i], yy = grid.y[j], r = grid.radius(i,j);
      // evaluate exact solution (though velocity info discarded for shelf)
      ierr = exactM(r, &alpha, &Drr, 1.0e-12, 0.0, 1); // CHKERRQ(ierr);
      if (ierr != 0) {
        verbPrintf(1,grid.com,
                   "exactM() not successful; gives alpha = %f (m/a), at\n"
                   "   i,j,xx,yy,r = %d,%d,%f,%f,%f\n",
                   alpha*secpera, i, j, xx, yy, r);
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
        vel_bc(i,j).u = myu; 
        vel_bc(i,j).v = myv;
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
        mask[i][j] = MASK_OCEAN_AT_TIME_0;
      }
    }
  }  
  ierr = vh.end_access(); CHKERRQ(ierr);    
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);    
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vel_bc.end_access(); CHKERRQ(ierr);

  ierr = stress_balance->set_boundary_conditions(vMask, vel_bc); CHKERRQ(ierr);

  // communicate what we have set that has ghosts
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);
  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);
  ierr = vbed.beginGhostComm(); CHKERRQ(ierr);
  ierr = vbed.endGhostComm(); CHKERRQ(ierr);
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactSSAModel::reportErrors() {
  PetscErrorCode  ierr;
  PetscScalar exactmaxu, maxvecerr = 0.0, avvecerr = 0.0, 
              avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  PetscScalar gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0,
              gmaxuerr = 0.0, gmaxverr = 0.0;

  if (config.get_flag("do_pseudo_plastic_till")) {
    ierr = verbPrintf(1,grid.com, 
          "WARNING: numerical errors not valid for pseudo-plastic till\n"); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com, 
          "NUMERICAL ERRORS in velocity relative to exact solution:\n"); CHKERRQ(ierr);

  IceModelVec2V *vel_ssa;
  ierr = stress_balance->get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr);

  ierr = vel_ssa->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2, uexact, vexact;
      PetscScalar myx = grid.x[i], myy = grid.y[j], myr = grid.radius(i,j);
      // eval exact solution
      if (test == 'I') {
        exactI(m_schoof, myx, myy, &junk1, &junk2, &uexact, &vexact); 
      } else if (test == 'J') {
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
      const PetscScalar uerr = PetscAbsReal((*vel_ssa)(i,j).u - uexact);
      const PetscScalar verr = PetscAbsReal((*vel_ssa)(i,j).v - vexact);
      avuerr = avuerr + uerr;      
      avverr = avverr + verr;      
      maxuerr = PetscMax(maxuerr,uerr);
      maxverr = PetscMax(maxverr,verr);
      const PetscScalar vecerr = sqrt(uerr * uerr + verr * verr);
      maxvecerr = PetscMax(maxvecerr,vecerr);
      avvecerr = avvecerr + vecerr;
    }
  }
  ierr = vel_ssa->end_access(); CHKERRQ(ierr);
   
  ierr = PetscGlobalMax(&maxuerr, &gmaxuerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxverr, &gmaxverr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avuerr, &gavuerr, grid.com); CHKERRQ(ierr);
  gavuerr = gavuerr/(grid.Mx*grid.My);
  ierr = PetscGlobalSum(&avverr, &gavverr, grid.com); CHKERRQ(ierr);
  gavverr = gavverr/(grid.Mx*grid.My);
  ierr = PetscGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.Mx*grid.My);

  // NetCDF error report code:
  int start;
  char filename[TEMPORARY_STRING_LENGTH];
  PetscTruth netcdf_report;
  NCTimeseries err;
  ierr = PetscOptionsGetString(PETSC_NULL, "-report_file", filename,
			       TEMPORARY_STRING_LENGTH, &netcdf_report); CHKERRQ(ierr);

  if (netcdf_report) {
    ierr = verbPrintf(2,grid.com, "Also writing errors to '%s'...\n", filename);
    CHKERRQ(ierr);

    // Find the number of records in this file:
    NCTool nc(grid.com, grid.rank);
    ierr = nc.open_for_writing(filename); CHKERRQ(ierr);
    ierr = nc.get_dim_length("N", &start); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = global_attributes.write(filename); CHKERRQ(ierr);

    // Write the dimension variable:
    err.init("N", "N", grid.com, grid.rank);
    ierr = err.write(filename, (size_t)start, (double)(start + 1), NC_INT); CHKERRQ(ierr);

    // Always write grid parameters:
    err.short_name = "dx";
    ierr = err.set_units("meters"); CHKERRQ(ierr);
    ierr = err.write(filename, (size_t)start, grid.dx); CHKERRQ(ierr);
    err.short_name = "dy";
    ierr = err.write(filename, (size_t)start, grid.dy); CHKERRQ(ierr);
    err.short_name = "dz";
    ierr = err.write(filename, (size_t)start, grid.dzMAX); CHKERRQ(ierr);
    err.short_name = "dzb";
    ierr = err.write(filename, (size_t)start, grid.dzbMAX); CHKERRQ(ierr);

    // Always write the test name:
    err.reset();
    err.short_name = "test";
    ierr = err.write(filename, (size_t)start, (double)test, NC_BYTE); CHKERRQ(ierr);
  }

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

    if (netcdf_report) {
      err.reset();
      err.short_name = "maximum_velocity";
      ierr = err.set_units("meters/second"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("meters/year"); CHKERRQ(ierr);
      ierr = err.write(filename, (size_t)start, gmaxvecerr); CHKERRQ(ierr);

      err.short_name = "maximum_u";
      ierr = err.write(filename, (size_t)start, gmaxuerr); CHKERRQ(ierr);
      err.short_name = "average_u";
      ierr = err.write(filename, (size_t)start, gavuerr); CHKERRQ(ierr);

      err.reset();
      err.short_name = "relative_velocity";
      ierr = err.set_units("percent"); CHKERRQ(ierr);
      ierr = err.write(filename, (size_t)start, (gavvecerr/exactmaxu)*100.0); CHKERRQ(ierr);
    }

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

    if (netcdf_report) {
      err.reset();
      err.short_name = "maximum_velocity";
      ierr = err.set_units("meters/second"); CHKERRQ(ierr);
      ierr = err.set_glaciological_units("meters/year"); CHKERRQ(ierr);
      ierr = err.write(filename, (size_t)start, gmaxvecerr); CHKERRQ(ierr);

      err.short_name = "maximum_u";
      ierr = err.write(filename, (size_t)start, gmaxuerr); CHKERRQ(ierr);
      err.short_name = "average_u";
      ierr = err.write(filename, (size_t)start, gavuerr); CHKERRQ(ierr);

      err.short_name = "maximum_v";
      ierr = err.write(filename, (size_t)start, gmaxverr); CHKERRQ(ierr);
      err.short_name = "average_v";
      ierr = err.write(filename, (size_t)start, gavverr); CHKERRQ(ierr);

      err.reset();
      err.short_name = "relative_velocity";
      ierr = err.set_units("percent"); CHKERRQ(ierr);
      ierr = err.write(filename, (size_t)start, (gavvecerr/exactmaxu)*100.0); CHKERRQ(ierr);
    }
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
  IceModelVec2V *vel_2d;
  IceModelVec3 *u3, *v3, *w3;

  ierr = stress_balance->get_advective_2d_velocity(vel_2d); CHKERRQ(ierr);
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr); 

  ierr = vel_2d->begin_access(); CHKERRQ(ierr);
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2;
      PetscScalar myx = grid.x[i], myy = grid.y[j],
        myr = grid.radius(i,j);
      if (test == 'I') {
        exactI(m_schoof, myx, myy, &junk1, &junk2, &(*vel_2d)(i,j).u, &(*vel_2d)(i,j).v); 
      } else if (test == 'J') {
        ierr = exactJ(myx, myy, &junk1, &junk2, &(*vel_2d)(i,j).u, &(*vel_2d)(i,j).v); CHKERRQ(ierr);
      } else if (test == 'M') {
        PetscScalar alpha;
        ierr = exactM(myr, &alpha, &junk1, 1.0e-12, 0.0, 1); //CHKERRQ(ierr);
        if (ierr != 0) {
          verbPrintf(1,grid.com,"exactM() not successful; gives alpha = %f (m/a), at\n",
                     alpha*secpera);
          verbPrintf(1,grid.com,"   i,j,xx,yy,r = %d,%d,%f,%f,%f\n", i,j,myx,myy,myr);
        }
        if (myr > 1.0) {
          (*vel_2d)(i,j).u = alpha * (myx / myr);
          (*vel_2d)(i,j).v = alpha * (myy / myr);
        } else {
          (*vel_2d)(i,j).u = 0.0;
          (*vel_2d)(i,j).v = 0.0;
        }
      } else {
        SETERRQ(1,"only tests I,J,M supported in IceExactSSAModel");
      }
      // now fill in velocity at depth
      for (PetscInt k=0; k<grid.Mz; k++) {
        ierr = u3->setColumn(i,j,(*vel_2d)(i,j).u); CHKERRQ(ierr);
        ierr = v3->setColumn(i,j,(*vel_2d)(i,j).v); CHKERRQ(ierr);
        // leave w alone for now; will be filled by call to broadcastSSAVelocity()
      }
    }
  }
  ierr = vel_2d->end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::run() {
  PetscErrorCode  ierr;

  ierr = verbPrintf(2,grid.com, "running Test %c ...\n", test); CHKERRQ(ierr);

  // does user want to turn off actual numerical evolution (and save exact solution)?
  bool flag;
  ierr = PISMOptionsIsSet("-eo", flag); CHKERRQ(ierr);
  if (flag) {
    exactOnly = PETSC_TRUE;
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
    ierr = stress_balance->update(false); CHKERRQ(ierr);
  }

  // report on result of computation (i.e. to standard out and to viewers)
  ierr = computeMax3DVelocities(); CHKERRQ(ierr); 
  ierr = summary(true); CHKERRQ(ierr);
  ierr = update_viewers(); CHKERRQ(ierr);

  return 0;
}
