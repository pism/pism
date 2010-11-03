// Copyright (C) 2004--2010 Constantine Khroulev, Ed Bueler and Jed Brown
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

#include "SSAFD.hh"

//! \brief Initialize the SSA solver.
PetscErrorCode SSAFD::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = ShallowStressBalance::init(vars); CHKERRQ(ierr);

  ierr = velocity.set(0.0); CHKERRQ(ierr);

  ierr = allocate(); CHKERRQ(ierr); 

  return 0;
}


//! \brief Solves the SSA if fast == false.
PetscErrorCode SSAFD::update(bool fast) {
  PetscErrorCode ierr;

  if (fast)
    return 0;

  ierr = compute_hardav_staggered(); CHKERRQ(ierr);
  
  ierr = assemble_rhs(); CHKERRQ(ierr);

  ierr = solve(); CHKERRQ(ierr);

  return 0;
}

//! \brief Allocate objects specific to the SSAFD object.
PetscErrorCode SSAFD::allocate() {
  PetscErrorCode ierr;

  // mimic IceGrid::createDA() with TRANSPOSE :
  PetscInt dof=2, stencil_width=1;
  ierr = DACreate2d(grid.com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    grid.My, grid.Mx,
                    grid.Ny, grid.Nx,
                    dof, stencil_width,
                    grid.procs_y, grid.procs_x,
                    &SSADA); CHKERRQ(ierr);

  ierr = DACreateGlobalVector(SSADA, &SSAX); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &SSARHS); CHKERRQ(ierr);

  ierr = DAGetMatrix(SSADA, MATMPIAIJ, &SSAStiffnessMatrix); CHKERRQ(ierr);

  ierr = KSPCreate(grid.com, &SSAKSP); CHKERRQ(ierr);
  // the default PC type somehow is ILU, which now fails (?) while block jacobi
  //   seems to work; runtime options can override (see test J in vfnow.py)
  PC pc;
  ierr = KSPGetPC(SSAKSP,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(SSAKSP); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode SSAFD::deallocate() {
  PetscErrorCode ierr;

  ierr = KSPDestroy(SSAKSP); CHKERRQ(ierr);
  ierr = MatDestroy(SSAStiffnessMatrix); CHKERRQ(ierr);
  ierr = VecDestroy(SSAX); CHKERRQ(ierr);
  ierr = VecDestroy(SSARHS); CHKERRQ(ierr);
  ierr = DADestroy(SSADA);CHKERRQ(ierr);
  
  return 0;
}

//! \brief Update the SSA solution.
PetscErrorCode SSAFD::update(bool fast) {
  PetscErrorCode ierr;

  if (fast)
    return 0;

  ierr = solve(); CHKERRQ(ierr); 

  ierr = compute_basal_frictional_heating(basal_frictional_heating); CHKERRQ(ierr);

  ierr = compute_D2(D2); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the D2 term (for the strain heating computation).
PetscErrorCode SSAFD::compute_D2(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = mask.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (mask->value(i,j) == MASK_DRAGGING_SHEET) {
        const PetscScalar 
          u_x   = (velocity(i+1,j).u - velocity(i-1,j).u)/(2*dx),
          u_y   = (velocity(i,j+1).u - velocity(i,j-1).u)/(2*dy),
          v_x   = (velocity(i+1,j).v - velocity(i-1,j).v)/(2*dx),
          v_y   = (velocity(i,j+1).v - velocity(i,j-1).v)/(2*dy),
          D2ssa = PetscSqr(u_x) + PetscSqr(v_y) + u_x * v_y
          + PetscSqr(0.5*(u_y + v_x));

        result(i,j) = D2ssa;
      } else {
        result(i,j) = 0.0;
      }
    }
  }

  ierr = mask.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the basal frictional heating.
PetscErrorCode SSAFD::compute_basal_frictional_heating(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (mask->is_floating(i,j)) {
        result(i,j) = 0.0;
      } else {
        const PetscScalar 
          C = basal->drag((*tauc)(i,j), velocity(i,j).u, velocity(i,j).v),
	  basal_stress_x = - C * velocity(i,j).u,
	  basal_stress_y = - C * velocity(i,j).v;
        result = basal_stress_x * velocity(i,j).u - basal_stress_y * velocity(i,j).v;
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);

  return 0;
}


//! \brief Computes vertically-averaged ice hardness on the staggered grid.
PetscErrorCode SSAFD::compute_hardav_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  PetscScalar *E, *E_ij, *E_offset;

  E = new PetscScalar[grid.Mz];

  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = enthalpy->getInternalColumn(i,j,&E_ij); CHKERRQ(ierr);
      for (PetscInt o=0; o<2; o++) {
        const PetscInt oi = 1-o, oj=o;  
        const PetscScalar H = 0.5 * ((*thickness)(i,j) + (*thickness)(i+oi,j+oj));

        if (H == 0) {
          result(i,j,o) = -1e6; // an obviously impossible value
          continue;
        }

        ierr = enthalpy->getInternalColumn(i+oi,j+oj,&E_offset); CHKERRQ(ierr);
        // build a column of enthalpy values a the current location:
        for (int k = 0; k < grid.Mz; ++k) {
          E[k] = 0.5 * (E_ij[k] + E_offset[k]);
        }
        
        result(i,j,o) = ice->averagedHardness_from_enth(H, grid.kBelowHeight(H),
                                                        grid.zlevels, E); CHKERRQ(ierr); 
      } // o
    }   // j
  }     // i

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  delete [] E;

  return 0;
}


//! \brief Assemble the right-hand-side vector.
PetscErrorCode SSAFD::assemble_rhs(Vec rhs) {
  PetscErrorCode ierr;
  PetscScalar     **taudx, **taudy;
  PISMVector2     **rhs_uv;

  // next constant not too sensitive, but must match value in assembleSSAMatrix():
  const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream;

  ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

  // get driving stress components
  ierr = compute_driving_stress(vWork2d[0],vWork2d[1]); CHKERRQ(ierr);

  ierr = vWork2d[0].get_array(taudx); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(taudy); CHKERRQ(ierr);
  ierr = bc_locations->begin_access(); CHKERRQ(ierr);
  ierr = DAVecGetArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);

  if (set_bc) {
    ierr = vel_bc->begin_access(); CHKERRQ(ierr);
  }

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (set_bc && (bc_locations->value(i,j) == MASK_BC)) {
        rhs_uv[i][j].u = scaling * vel_bc(i,j).u;
        rhs_uv[i][j].v = scaling * vel_bc(i,j).v;
      } else {
	// usual case: use already computed driving stress
        rhs_uv[i][j].u = taudx[i][j];
        rhs_uv[i][j].v = taudy[i][j];
      }
    }
  }

  if (set_bc) {
    ierr = vel_bc->end_access(); CHKERRQ(ierr);
  }

  ierr = DAVecRestoreArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);
  ierr = bc_locations->end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);

  return 0;
}

//! \brief Assemble the SSA matrix.
PetscErrorCode SSAFD::assemble_matrix(bool include_basal_shear, Mat A) {
  PetscErrorCode ierr;
  
  return 0;
}

//! \brief Compute the norm of nuH and the change in nuH.
PetscErrorCode SSAFD::compute_nuH_norm(IceModelVec2Stag nuH,
                                       IceModelVec2Stag nuH_old,
                                       PetscReal &norm,
                                       PetscReal &norm_change) {
  PetscErrorCode ierr;

  PetscReal nuNorm[2], nuChange[2];

  const PetscScalar area = grid.dx * grid.dy;
#define MY_NORM     NORM_1

  // Test for change in nu
  ierr = nuH_old.add(-1, nuH); CHKERRQ(ierr);

  ierr = nuH_old.norm_all(MY_NORM, nuChange[0], nuChange[1]); CHKERRQ(ierr);
  ierr =     nuH.norm_all(MY_NORM, nuNorm[0],   nuNorm[1]);   CHKERRQ(ierr);

  nuChange[0] *= area;
  nuChange[1] *= area;
  nuNorm[0] *= area;
  nuNorm[1] *= area;

  normChange = sqrt(PetscSqr(nuChange[0]) + PetscSqr(nuChange[1]));
  norm = sqrt(PetscSqr(nuNorm[0]) + PetscSqr(nuNorm[1]));
  
  return 0;
}

//! \brief Compute the product of ice thickness and effective viscosity.
PetscErrorCode SSAFD::compute_nuH_staggered(IceModelVec2Stag &hardness,
                                            IceModelVec2Stag &result,
                                            PetscReal epsilon) {
  PetscErrorCode ierr;
  PISMVector2 **uv;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = velocity.get_array(uv); CHKERRQ(ierr);
  ierr = hardness.begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  const PetscScalar dx = grid.dx, dy = grid.dy;

  for (PetscInt o=0; o<2; ++o) {
    const PetscInt oi = 1 - o, oj=o;
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

        const PetscScalar H = 0.5 * ((*thickness)(i,j) + (*thickness)(i+oi,j+oj));

        if (H < ssaStrengthExtend.min_thickness_for_extension()) {
          // Extends strength of SSA (i.e. nuH coeff) into the ice free region.
          //  Does not add or subtract ice mass.
          result(i,j,o) = ssaStrengthExtend.notional_strength();
        } else {
          PetscScalar u_x, u_y, v_x, v_y;
          // Check the offset to determine how to differentiate velocity
          if (o == 0) {
            u_x = (uv[i+1][j].u - uv[i][j].u) / dx;
            u_y = (uv[i][j+1].u + uv[i+1][j+1].u - uv[i][j-1].u - uv[i+1][j-1].u) / (4*dy);
            v_x = (uv[i+1][j].v - uv[i][j].v) / dx;
            v_y = (uv[i][j+1].v + uv[i+1][j+1].v - uv[i][j-1].v - uv[i+1][j-1].v) / (4*dy);
          } else {
            u_x = (uv[i+1][j].u + uv[i+1][j+1].u - uv[i-1][j].u - uv[i-1][j+1].u) / (4*dx);
            u_y = (uv[i][j+1].u - uv[i][j].u) / dy;
            v_x = (uv[i+1][j].v + uv[i+1][j+1].v - uv[i-1][j].v - uv[i-1][j+1].v) / (4*dx);
            v_y = (uv[i][j+1].v - uv[i][j].v) / dy;
          }

          const PetscScalar hardav = hardness(i,j,o);
          result(i,j,o) = H * ice->effectiveViscosity(hardav, u_x, u_y, v_x, v_y);

          if (! finite(result(i,j,o)) || false) {
            ierr = PetscPrintf(grid.com, "nuH[%d][%d][%d] = %e\n", o, i, j, result(i,j,o));
              CHKERRQ(ierr); 
            ierr = PetscPrintf(grid.com, "  u_x, u_y, v_x, v_y = %e, %e, %e, %e\n", 
                               u_x, u_y, v_x, v_y);
              CHKERRQ(ierr);
          }
          
          // We ensure that nuH is bounded below by a positive constant.
          result(i,j,o) += epsilon;
        } // end of if (H < ssaStrengthExtend.min_thickness_for_extension()) { ... } else {
      } // j
    } // i
  } // o

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = hardness.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);

  // Some communication
  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);
  
  return 0;
}

//! \brief Compute the driving stress.
/*!
 * FIXME: this method should use the surface gradient computed "the right way"
 * instead of doing its own thing here.
 */
PetscErrorCode SSAFD::compute_driving_stress(IceModelVec2S &taudx,
                                             IceModelVec2S &taudy) {
  PetscErrorCode ierr;
  
  const PetscScalar n       = ice->exponent(), // frequently n = 3
                    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
                    invpow  = 1.0 / etapow,  // = 3/8
                    dinvpow = (- n - 2.0) / (2.0 * n + 2.0); // = -5/8
  const PetscScalar minThickEtaTransform = 5.0; // m
  const PetscScalar dx=grid.dx, dy=grid.dy;

  bool compute_surf_grad_inward_ssa = config.get_flag("compute_surf_grad_inward_ssa");
  bool use_eta = (config.get_string("surface_gradient_method") == "eta");

  ierr =   surface->begin_access();    CHKERRQ(ierr);
  ierr = thickness->begin_access();  CHKERRQ(ierr);
  ierr =       bed->begin_access();  CHKERRQ(ierr);
  ierr =      mask->begin_access();  CHKERRQ(ierr);

  ierr = taudx.begin_access(); CHKERRQ(ierr);
  ierr = taudy.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscReal thk = (*thickness)(i,j);
      const PetscScalar pressure = ice->rho * standard_gravity * thk;
      if (pressure <= 0.0) {
        vtaudx(i,j) = 0.0;
        vtaudy(i,j) = 0.0;
      } else {
        PetscScalar h_x = 0.0, h_y = 0.0;
        // FIXME: we need to handle grid periodicity correctly.
        if (vMask.is_grounded(i,j) && (use_eta == true)) {
	        // in grounded case, differentiate eta = H^{8/3} by chain rule
          if (thk > 0.0) {
            const PetscScalar myH = (thk < minThickEtaTransform) ?
	                                  minThickEtaTransform : thk;
            const PetscScalar eta = pow(myH, etapow), factor = invpow * pow(eta, dinvpow);
            h_x = factor * (pow((*thickness)(i+1,j),etapow) - pow((*thickness)(i-1,j),etapow)) / (2*dx);
            h_y = factor * (pow((*thickness)(i,j+1),etapow) - pow((*thickness)(i,j-1),etapow)) / (2*dy);
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
        }

        taudx(i,j) = - pressure * h_x;
        taudy(i,j) = - pressure * h_y;
      }
    }
  }

  ierr =       bed->end_access(); CHKERRQ(ierr);
  ierr =   surface->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr =      mask->end_access(); CHKERRQ(ierr);
  ierr =      taudx.end_access(); CHKERRQ(ierr);
  ierr =      taudy.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFD::solve() {
  PetscErrorCode ierr;

  return 0;
}

