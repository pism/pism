/* Copyright (C) 2014 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "SIAFD_FEvoR.hh"

#include "PISMVars.hh"
#include "flowlaws.hh"
#include "PISMBedSmoother.hh"
#include "enthalpyConverter.hh"

namespace pism {

SIAFD_FEvoR::SIAFD_FEvoR(IceGrid &g, EnthalpyConverter &e, const Config &c)
  : SIAFD(g, e, c) {
  // empty
}

SIAFD_FEvoR::~SIAFD_FEvoR() {
  // empty
}

PetscErrorCode SIAFD_FEvoR::init(Vars &vars) {
  PetscErrorCode ierr;

  ierr = SIAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "  [using the enhancement factor computed using FEvoR]\n"); CHKERRQ(ierr);

  m_variables = &vars;

  return 0;
}

PetscErrorCode SIAFD_FEvoR::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  if (config.get_flag("siafd_fevor_use_constant_slope")) {
    double slope = (config.get("siafd_fevor_surface_slope_degrees") / 180.0) * M_PI;

    ierr = h_x.set(slope); CHKERRQ(ierr);
    ierr = h_y.set(0.0); CHKERRQ(ierr);
  } else {
    ierr = SIAFD::compute_surface_gradient(h_x, h_y); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode SIAFD_FEvoR::compute_diffusive_flux(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y,
                                                   IceModelVec2Stag &result, bool fast) {
  PetscErrorCode  ierr;
  IceModelVec2S &thk_smooth = work_2d[0],
    &theta = work_2d[1];

  bool full_update = (fast == false);

  ierr = result.set(0.0); CHKERRQ(ierr);

  std::vector<double> delta_ij(grid.Mz);

  double ice_grain_size = config.get("ice_grain_size");

  bool compute_grain_size_using_age = config.get_flag("compute_grain_size_using_age");

  // some flow laws use grain size, and even need age to update grain size
  if (compute_grain_size_using_age && (!config.get_flag("do_age"))) {
    PetscPrintf(grid.com,
                "PISM ERROR in SIAFD::compute_diffusive_flux(): do_age not set but\n"
                "age is needed for grain-size-based flow law ...  ENDING! ...\n\n");
    PISMEnd();
  }

  const bool use_age = (IceFlowLawUsesGrainSize(flow_law) &&
                        compute_grain_size_using_age &&
                        config.get_flag("do_age"));

  // get "theta" from Schoof (2003) bed smoothness calculation and the
  // thickness relative to the smoothed bed; each IceModelVec2S involved must
  // have stencil width WIDE_GHOSTS for this too work
  ierr = bed_smoother->get_theta(*surface, &theta); CHKERRQ(ierr);

  ierr = bed_smoother->get_smoothed_thk(*surface, *thickness, *mask,
                                        &thk_smooth); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(theta);
  list.add(thk_smooth);
  list.add(result);

  list.add(h_x);
  list.add(h_y);

  double *age_ij, *age_offset;
  if (use_age) {
    list.add(*age);
  }

  if (full_update) {
    list.add(delta[0]);
    list.add(delta[1]);
  }

  double *E_ij, *E_offset;
  list.add(*enthalpy);

  // new code
  IceModelVec3 *enhancement_factor = dynamic_cast<IceModelVec3*>(m_variables->get("enhancement_factor"));
  if (enhancement_factor == NULL) {
    SETERRQ(grid.com, 1, "enhancement_factor is not available");
  }

  double *EF_ij, *EF_offset;
  list.add(*enhancement_factor);
  // end of new code

  double my_D_max = 0.0;
  for (int o=0; o<2; o++) {
    for (PointsWithGhosts p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // staggered point: o=0 is i+1/2, o=1 is j+1/2, (i,j) and (i+oi,j+oj)
      //   are regular grid neighbors of a staggered point:
      const int oi = 1 - o, oj = o;

      const double
        thk = 0.5 * (thk_smooth(i,j) + thk_smooth(i+oi,j+oj));

      // zero thickness case:
      if (thk == 0.0) {
        result(i,j,o) = 0.0;
        if (full_update) {
          ierr = delta[o].setColumn(i, j, 0.0); CHKERRQ(ierr);
        }
        continue;
      }

      if (use_age) {
        ierr = age->getInternalColumn(i, j, &age_ij); CHKERRQ(ierr);
        ierr = age->getInternalColumn(i+oi, j+oj, &age_offset); CHKERRQ(ierr);
      }

      ierr = enthalpy->getInternalColumn(i, j, &E_ij); CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i+oi, j+oj, &E_offset); CHKERRQ(ierr);

      // new code
      ierr = enhancement_factor->getInternalColumn(i, j, &EF_ij); CHKERRQ(ierr);
      ierr = enhancement_factor->getInternalColumn(i+oi, j+oj, &EF_offset); CHKERRQ(ierr);
      // end of new code

      const double slope = (o==0) ? h_x(i,j,o) : h_y(i,j,o);
      const int      ks = grid.kBelowHeight(thk);
      const double   alpha =
        sqrt(PetscSqr(h_x(i,j,o)) + PetscSqr(h_y(i,j,o)));
      const double theta_local = 0.5 * (theta(i,j) + theta(i+oi,j+oj));

      double  Dfoffset = 0.0;  // diffusivity for deformational SIA flow
      for (int k = 0; k <= ks; ++k) {
        double depth = thk - grid.zlevels[k]; // FIXME issue #15
        // pressure added by the ice (i.e. pressure difference between the
        // current level and the top of the column)
        const double pressure = EC.getPressureFromDepth(depth);

        double flow;
        if (use_age) {
          ice_grain_size = grainSizeVostok(0.5 * (age_ij[k] + age_offset[k]));
        }
        // If the flow law does not use grain size, it will just ignore it,
        // no harm there
        double E = 0.5 * (E_ij[k] + E_offset[k]);
        flow = flow_law->flow(alpha * pressure, E, pressure, ice_grain_size);

        // new code
        // compute the enhancement factor on the staggered grid
        double EF = 0.5 * (EF_ij[k] + EF_offset[k]);

        delta_ij[k] = EF * theta_local * 2.0 * pressure * flow;
        // end of new code

        if (k > 0) { // trapezoidal rule
          const double dz = grid.zlevels[k] - grid.zlevels[k-1];
          Dfoffset += 0.5 * dz * ((depth + dz) * delta_ij[k-1] + depth * delta_ij[k]);
        }
      }
      // finish off D with (1/2) dz (0 + (H-z[ks])*delta_ij[ks]), but dz=H-z[ks]:
      const double dz = thk - grid.zlevels[ks];
      Dfoffset += 0.5 * dz * dz * delta_ij[ks];

      // Override diffusivity at the edges of the domain. (At these
      // locations PISM uses ghost cells *beyond* the boundary of
      // the computational domain. This does not matter if the ice
      // does not extend all the way to the domain boundary, as in
      // whole-ice-sheet simulations. In a regional setup, though,
      // this adjustment lets us avoid taking very small time-steps
      // because of the possible thickness and bed elevation
      // "discontinuities" at the boundary.)
      if (i < 0 || i >= grid.Mx-1 || j < 0 || j >= grid.My-1)
        Dfoffset = 0.0;

      my_D_max = PetscMax(my_D_max, Dfoffset);

      // vertically-averaged SIA-only flux, sans sliding; note
      //   result(i,j,0) is  u  at E (east)  staggered point (i+1/2,j)
      //   result(i,j,1) is  v  at N (north) staggered point (i,j+1/2)
      result(i,j,o) = - Dfoffset * slope;

      // if doing the full update, fill the delta column above the ice and
      // store it:
      if (full_update) {
        for (unsigned int k = ks + 1; k < grid.Mz; ++k) {
          delta_ij[k] = 0.0;
        }
        ierr = delta[o].setInternalColumn(i,j,&delta_ij[0]); CHKERRQ(ierr);
      }
    }
  } // o-loop

  ierr = GlobalMax(grid.com, &my_D_max, &D_max); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
