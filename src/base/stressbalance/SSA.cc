// Copyright (C) 2004--2012 Constantine Khroulev, Ed Bueler, Jed Brown, Torsten Albrecht
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

#include "SSA.hh"
#include "Mask.hh"
#include "basal_resistance.hh"
#include "PISMVars.hh"
#include "PISMProf.hh"
#include "pism_options.hh"
#include "flowlaw_factory.hh"
#include "PIO.hh"

SSA::SSA(IceGrid &g, IceBasalResistancePlasticLaw &b,
         EnthalpyConverter &e,
         const NCConfigVariable &c)
  : ShallowStressBalance(g, b, e, c)
{
  mask = NULL;
  thickness = NULL;
  tauc = NULL;
  surface = NULL;
  bed = NULL;
  enthalpy = NULL;
  driving_stress_x = NULL;
  driving_stress_y = NULL;

  strength_extension = new SSAStrengthExtension(config);
  allocate();
}


//! \brief Initialize a generic regular-grid SSA solver.
PetscErrorCode SSA::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = ShallowStressBalance::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"* Initializing the SSA stress balance...\n"); CHKERRQ(ierr);

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  tauc = dynamic_cast<IceModelVec2S*>(vars.get("tauc"));
  if (tauc == NULL) SETERRQ(grid.com, 1, "tauc is not available");

  surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  driving_stress_x = dynamic_cast<IceModelVec2S*>(vars.get("ssa_driving_stress_x"));
  driving_stress_y = dynamic_cast<IceModelVec2S*>(vars.get("ssa_driving_stress_y"));
  if( (driving_stress_x==NULL) && (driving_stress_y==NULL) ) {
    if(surface == NULL) {
      SETERRQ(grid.com, 1, "neither surface_altitude nor ssa_driving_stress_x/y is available");      
    }
  } else if(surface !=NULL){
    SETERRQ(grid.com, 1, "at most one of surface_altitude or ssa_driving_stress_x/y may be specified");    
  } else if( (driving_stress_x==NULL) || (driving_stress_y==NULL) ) {
    SETERRQ(grid.com, 1, "both of ssa_driving_stress_x/y must be specified if one is");
  }

  bed = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (bed == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");


  // Check if PISM is being initialized from an output file from a previous run
  // and read the initial guess (unless asked not to).
  bool i_set;
  string filename;
  ierr = PISMOptionsString("-i", "PISM input file",
                           filename, i_set); CHKERRQ(ierr);

  if (i_set) {
    bool dont_read_initial_guess, u_ssa_found, v_ssa_found;
    unsigned int start;
    PIO nc(grid.com, grid.rank, "netcdf3");

    ierr = PISMOptionsIsSet("-dontreadSSAvels", dont_read_initial_guess); CHKERRQ(ierr);

    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_var("u_ssa", u_ssa_found); CHKERRQ(ierr); 
    ierr = nc.inq_var("v_ssa", v_ssa_found); CHKERRQ(ierr); 
    ierr = nc.inq_nrecords(start); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr); 
    start -= 1;

    if (u_ssa_found && v_ssa_found &&
        (! dont_read_initial_guess)) {
      ierr = verbPrintf(3,grid.com,"Reading u_ssa and v_ssa...\n"); CHKERRQ(ierr);

      ierr = velocity.read(filename.c_str(), start); CHKERRQ(ierr); 
    }
    
  } else {
    ierr = velocity.set(0.0); CHKERRQ(ierr); // default initial guess
  }

  if (config.get_flag("ssa_dirichlet_bc")) {
    bc_locations = dynamic_cast<IceModelVec2Int*>(vars.get("bcflag"));
    if (bc_locations == NULL) SETERRQ(grid.com, 1, "bc_locations is not available");
    
    vel_bc = dynamic_cast<IceModelVec2V*>(vars.get("vel_ssa_bc"));
    if (vel_bc == NULL) SETERRQ(grid.com, 1, "vel_ssa_bc is not available");
  }

  event_ssa = grid.profiler->create("ssa_update", "time spent solving the SSA");

  return 0;
}

//! \brief Allocate objects which any SSA solver would use.
PetscErrorCode SSA::allocate() {
  PetscErrorCode ierr;

  ierr = taud.create(grid, "taud", false); CHKERRQ(ierr);
  ierr = taud.set_attrs("diagnostic",
                        "X-component of the driving shear stress at the base of ice",
                        "Pa", "", 0); CHKERRQ(ierr);
  ierr = taud.set_attrs("diagnostic",
                        "Y-component of the driving shear stress at the base of ice",
                        "Pa", "", 1); CHKERRQ(ierr);

  ierr = velocity_old.create(grid, "velocity_old", true); CHKERRQ(ierr);
  ierr = velocity_old.set_attrs("internal",
                                "old SSA velocity field; used for re-trying with a different epsilon",
                                "m s-1", ""); CHKERRQ(ierr);

  // override velocity metadata
  vector<string> long_names;
  long_names.push_back("SSA model ice velocity in the X direction");
  long_names.push_back("SSA model ice velocity in the Y direction");
  ierr = velocity.rename("_ssa",long_names,""); CHKERRQ(ierr);

  // mimic IceGrid::createDA() with TRANSPOSE :
  PetscInt dof=2, stencil_width=1;
  ierr = DMDACreate2d(grid.com,
                      DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      grid.My, grid.Mx,
                      grid.Ny, grid.Nx,
                      dof, stencil_width,
                      &grid.procs_y[0], &grid.procs_x[0],
                      &SSADA); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(SSADA, &SSAX); CHKERRQ(ierr);

  {
    IceFlowLawFactory ice_factory(grid.com, "ssa_", config, &EC);
    ice_factory.removeType(ICE_GOLDSBY_KOHLSTEDT);

    ierr = ice_factory.setType(config.get_string("ssa_flow_law").c_str()); CHKERRQ(ierr);

    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ierr = ice_factory.create(&flow_law); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSA::deallocate() {
  PetscErrorCode ierr;

  if (SSAX != PETSC_NULL) {
    ierr = VecDestroy(&SSAX); CHKERRQ(ierr);
  }

  if (SSADA != PETSC_NULL) {
    ierr = DMDestroy(&SSADA);CHKERRQ(ierr);
  }

  if (flow_law != NULL) {
    delete flow_law;
    flow_law = NULL;
  }

  return 0;
}


//! \brief Update the SSA solution.
PetscErrorCode SSA::update(bool fast) {
  PetscErrorCode ierr;

  if (fast)
    return 0;

  grid.profiler->begin(event_ssa);

  ierr = solve(); CHKERRQ(ierr); 

  ierr = compute_basal_frictional_heating(basal_frictional_heating); CHKERRQ(ierr);
  ierr = compute_D2(D2); CHKERRQ(ierr);

  ierr = compute_maximum_velocity(); CHKERRQ(ierr);

  grid.profiler->end(event_ssa);

  return 0;
}


//! \brief Compute the D2 term for the strain heating computation.
/*!
Documented in [\ref BBssasliding].
 */
PetscErrorCode SSA::compute_D2(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal dx = grid.dx, dy = grid.dy;

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      const PetscScalar 
          u_x   = (velocity(i+1,j).u - velocity(i-1,j).u)/(2*dx),
          u_y   = (velocity(i,j+1).u - velocity(i,j-1).u)/(2*dy),
          v_x   = (velocity(i+1,j).v - velocity(i-1,j).v)/(2*dx),
          v_y   = (velocity(i,j+1).v - velocity(i,j-1).v)/(2*dy);
      result(i,j) = PetscSqr(u_x) + PetscSqr(v_y) + u_x * v_y 
                      + PetscSqr(0.5*(u_y + v_x));
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);
  return 0;
}


//! \brief Compute eigenvalues of the horizontal, vertically-integrated strain rate tensor.
/*!
Calculates all components \f$D_{xx}, D_{yy}, D_{xy}=D_{yx}\f$ of the
vertically-averaged strain rate tensor \f$D\f$ [\ref SchoofStream].  Then computes
the eigenvalues \c result_e1 = (maximum eigenvalue), \c result_e1 = (minimum
eigenvalue).  Uses the provided thickness to make decisions (PIK) about computing
strain rates near calving front.

Though there are two eigenvalues, such do not form a vector, so the output is not
an IceModelVec2V, though it could be a std::vector<IceModelVec2S> or such.

Note that \c result_e1 >= \c result_e2, but there is no necessary relation between 
the magnitudes, and either principal strain rate could be negative or positive.

Result can be used in a calving law, for example in eigencalving (PIK).

FIXME:  need to answer: strain rates will be derived from SSA velocities. Is there ghost communication needed?
*/
PetscErrorCode SSA::compute_principal_strain_rates(IceModelVec2S &result_e1,
                                                   IceModelVec2S &result_e2) {
  PetscErrorCode ierr;
  PetscScalar    dx = grid.dx, dy = grid.dy;

  IceModelVec2S H = *thickness; // an alias
  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = H.begin_access(); CHKERRQ(ierr);
  ierr = result_e1.begin_access(); CHKERRQ(ierr);
  ierr = result_e2.begin_access(); CHKERRQ(ierr);

  Mask M;
  ierr = mask->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if ( velocity(i,j).u == 0.0 || velocity(i,j).v == 0.0) {
        result_e1(i,j) = 0.0;
        result_e2(i,j) = 0.0;
        continue;
      }

      const PetscInt
        M_e = mask->as_int(i + 1,j),
        M_w = mask->as_int(i - 1,j),
        M_n = mask->as_int(i,j + 1),
        M_s = mask->as_int(i,j - 1);

      //centered difference scheme; strain in units s-1
      PetscScalar
        u_x = (velocity(i+1,j).u - velocity(i-1,j).u) / (2 * dx),
        u_y = (velocity(i,j+1).u - velocity(i,j-1).u) / (2 * dy),
        v_x = (velocity(i+1,j).v - velocity(i-1,j).v) / (2 * dx),
        v_y = (velocity(i,j+1).v - velocity(i,j-1).v) / (2 * dy);



      //inward scheme at the ice-shelf
      //SSA velocity exists depending on mask (newly filled grid cells are not taken into account)

      if (M.ice_free(M_e)) {
        u_x = (velocity(i,j).u - velocity(i-1,j).u) / dx;
        v_x = (velocity(i,j).v - velocity(i-1,j).v) / dx;
      }
      if (M.ice_free(M_w)) {
        u_x = (velocity(i+1,j).u - velocity(i,j).u) / dx;
        v_x = (velocity(i+1,j).v - velocity(i,j).v) / dx;
      }
      if (M.ice_free(M_n)) {
        u_y = (velocity(i,j).u - velocity(i,j-1).u) / dy;
        v_y = (velocity(i,j).v - velocity(i,j-1).v) / dy;
      }
      if (M.ice_free(M_s)) {
        u_y = (velocity(i,j+1).u - velocity(i,j).u) / dy;
        v_y = (velocity(i,j+1).v - velocity(i,j).v) / dy;
      }

      // ice nose case
      if (M.ice_free(M_s) && M.ice_free(M_n)) {
        u_y = 0.0;
        v_y = 0.0;
      }
      if (M.ice_free(M_e) && M.ice_free(M_w)) {
        u_x = 0.0;
        v_x = 0.0;
      }

      const PetscScalar A = 0.5 * (u_x + v_y),  // A = (1/2) trace(D)
        B   = 0.5 * (u_x - v_y),
        Dxy = 0.5 * (v_x + u_y),  // B^2 = A^2 - u_x v_y
        q   = sqrt(PetscSqr(B) + PetscSqr(Dxy));
      result_e1(i,j) = A + q;
      result_e2(i,j) = A - q; // q >= 0 so e1 >= e2

    } // j
  }   // i

  ierr = velocity.end_access(); CHKERRQ(ierr);
  ierr = H.end_access(); CHKERRQ(ierr);
  ierr = result_e1.end_access(); CHKERRQ(ierr);
  ierr = result_e2.end_access(); CHKERRQ(ierr);

  ierr = mask->end_access(); CHKERRQ(ierr);
  return 0;
}


//! \brief Compute the basal frictional heating.
/*!
  Ice shelves have zero basal friction heating.
 */
PetscErrorCode SSA::compute_basal_frictional_heating(IceModelVec2S &result) {
  PetscErrorCode ierr;

  MaskQuery m(*mask);

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (m.ocean(i,j)) {
        result(i,j) = 0.0;
      } else {
        const PetscScalar 
          C = basal.drag((*tauc)(i,j), velocity(i,j).u, velocity(i,j).v),
	  basal_stress_x = - C * velocity(i,j).u,
	  basal_stress_y = - C * velocity(i,j).v;
        result(i,j) = - basal_stress_x * velocity(i,j).u - basal_stress_y * velocity(i,j).v;
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);
  return 0;
}

//! \brief Compute the driving stress.
/*!
Computes the driving stress at the base of the ice:
\f[ \tau_d = - \rho g H \nabla h \f]

If configuration parameter \c surface_gradient_method = \c eta then the surface
gradient \f$\nabla h\f$ is computed by the gradient of the transformed variable
\f$\eta= H^{(2n+2)/n}\f$ (frequently, \f$\eta= H^{8/3}\f$). The idea is that
this quantity is more regular at ice sheet margins, and so we get a better
surface gradient. When the thickness at a grid point is very small (below \c
minThickEtaTransform in the procedure), the formula is slightly modified to
give a lower driving stress. The transformation is not used in floating ice.
 */
PetscErrorCode SSA::compute_driving_stress(IceModelVec2V &result) {
  PetscErrorCode ierr;

  IceModelVec2S &thk = *thickness; // to improve readability (below)

  const PetscScalar n = flow_law->exponent(), // frequently n = 3
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,  // = 3/8
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0); // = -5/8
  const PetscScalar minThickEtaTransform = 5.0; // m
  const PetscScalar dx=grid.dx, dy=grid.dy;

  bool compute_surf_grad_inward_ssa = config.get_flag("compute_surf_grad_inward_ssa");
  PetscReal standard_gravity = config.get("standard_gravity"),
    ice_rho = config.get("ice_density");
  bool use_eta = (config.get_string("surface_gradient_method") == "eta");

  ierr =   surface->begin_access();    CHKERRQ(ierr);
  ierr =       bed->begin_access();  CHKERRQ(ierr);
  ierr =      mask->begin_access();  CHKERRQ(ierr);
  ierr =        thk.begin_access();  CHKERRQ(ierr);

  MaskQuery m(*mask);

  ierr = result.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar pressure = ice_rho * standard_gravity * thk(i,j); // FIXME issue #15
      if (pressure <= 0.0) {
        result(i,j).u = 0.0;
        result(i,j).v = 0.0;
      } else {
        PetscScalar h_x = 0.0, h_y = 0.0;
        // FIXME: we need to handle grid periodicity correctly.
        if (m.grounded(i,j) && (use_eta == true)) {
	        // in grounded case, differentiate eta = H^{8/3} by chain rule
          if (thk(i,j) > 0.0) {
            const PetscScalar myH = (thk(i,j) < minThickEtaTransform ?
                                     minThickEtaTransform : thk(i,j));
            const PetscScalar eta = pow(myH, etapow), factor = invpow * pow(eta, dinvpow);
            h_x = factor * (pow(thk(i+1,j),etapow) - pow(thk(i-1,j),etapow)) / (2*dx);
            h_y = factor * (pow(thk(i,j+1),etapow) - pow(thk(i,j-1),etapow)) / (2*dy);
          }
          // now add bed slope to get actual h_x,h_y
          // FIXME: there is no reason to assume user's bed is periodized
          h_x += bed->diff_x(i,j);
          h_y += bed->diff_y(i,j);
        } else {  // floating or eta transformation is not used
          if (compute_surf_grad_inward_ssa) {
            h_x = surface->diff_x_p(i,j);
            h_y = surface->diff_y_p(i,j);
          } else {
            h_x = surface->diff_x(i,j);
            h_y = surface->diff_y(i,j);
          }

	  // for floating shear margin we calculate inward scheme along ice free bedrock
	  bool ShearMarginE = (thk(i+1,j)<1.0 && (*bed)(i+1,j)>0.0),
	       ShearMarginW = (thk(i-1,j)<1.0 && (*bed)(i-1,j)>0.0),
	       ShearMarginN = (thk(i,j+1)<1.0 && (*bed)(i,j+1)>0.0),
	       ShearMarginS = (thk(i,j-1)<1.0 && (*bed)(i,j-1)>0.0);
	
	  bool shearMargin = (ShearMarginE || ShearMarginW || ShearMarginN || ShearMarginS);
	
	
	  if (shearMargin) {	
		
	    if (ShearMarginE && !ShearMarginW)
	      h_x = ((*surface)(i,j) - (*surface)(i-1,j)) / grid.dx;
	    else if (ShearMarginW && !ShearMarginE)
	      h_x = ((*surface)(i+1,j) - (*surface)(i,j)) / grid.dx;
	    else if (ShearMarginW && ShearMarginE)
	      h_x = 0.0;
		
	    if (ShearMarginN && !ShearMarginS)
	      h_y = ((*surface)(i,j) - (*surface)(i,j-1)) / grid.dy;
	    else if (ShearMarginS && !ShearMarginN)
	      h_y = ((*surface)(i,j+1) - (*surface)(i,j)) / grid.dy;
	    else if (ShearMarginN && ShearMarginS)
	      h_y = 0.0;
	  }
        }

        result(i,j).u = - pressure * h_x;
        result(i,j).v = - pressure * h_y;
      }
    }
  }

  ierr =        thk.end_access(); CHKERRQ(ierr);
  ierr =       bed->end_access(); CHKERRQ(ierr);
  ierr =   surface->end_access(); CHKERRQ(ierr);
  ierr =      mask->end_access(); CHKERRQ(ierr);
  ierr =     result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! \brief Compute maximum ice velocity; uses the mask to ignore values in
//!  ice-free areas.
PetscErrorCode SSA::compute_maximum_velocity() {
  PetscErrorCode ierr;
  PetscReal my_max_u = 0.0,
    my_max_v = 0.0;

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  
  MaskQuery m(*mask);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (m.icy(i, j)) {
        my_max_u = PetscMax(my_max_u, PetscAbs(velocity(i,j).u));
        my_max_v = PetscMax(my_max_v, PetscAbs(velocity(i,j).v));
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&my_max_u, &max_u, grid.com); CHKERRQ(ierr); 
  ierr = PISMGlobalMax(&my_max_v, &max_v, grid.com); CHKERRQ(ierr); 
  return 0;
}

PetscErrorCode SSA::stdout_report(string &result) {
  result = stdout_ssa;
  return 0;
}


//! \brief Set the initial guess of the SSA velocity.
PetscErrorCode SSA::set_initial_guess(IceModelVec2V &guess) {
  PetscErrorCode ierr;
  ierr = velocity.copy_from(guess); CHKERRQ(ierr);
  return 0;
}


void SSA::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["vel_ssa"] = velocity.get_metadata();
}


PetscErrorCode SSA::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "vel_ssa")) {
    ierr = velocity.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSA::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "vel_ssa")) {
    ierr = velocity.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

