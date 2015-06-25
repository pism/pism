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

#include <iostream>

#include "base/stressbalance/sia/SIAFD_FEvoR.hh"

#include "SIAFD.hh"
#include "PISMBedSmoother.hh"
#include "base/enthalpyConverter.hh"
#include "base/rheology/flowlaw_factory.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "base/util/Profiling.hh"


namespace pism {

namespace stressbalance {
  SIAFD_FEvoR::SIAFD_FEvoR(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e )
    : SIAFD(g, e) {

  // empty
}

SIAFD_FEvoR::~SIAFD_FEvoR() {
  // empty
}

void SIAFD_FEvoR::init() {

  SIAFD::init(); 

  m_log->message(2,
                    "  [using the enhancement factor computed using FEvoR]\n"); 

}

void SIAFD_FEvoR::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {

  if (m_config->get_boolean("sia_fevor_use_constant_slope")) {
    double slope = (m_config->get_double("sia_fevor_bed_slope_degrees") / 180.0) * M_PI;

    // We compute the surface slope using the fact that we are
    // modeling grounded ice, so surface = bed + thickness and
    // surface' = bed' + thickness'.

    const double dx = m_grid->dx(), dy = m_grid->dy();  // convenience

  const IceModelVec2S
    //    &h = *m_grid->variables().get_2d_scalar("surface_altitude"),
    &H = *m_grid->variables().get_2d_scalar("land_ice_thickness");


    IceModelVec::AccessList list;
    list.add(h_x);
    list.add(h_y);
    list.add(H);

    // h_x and h_y have to have ghosts
    assert(h_x.get_stencil_width() >= 1);
    assert(h_y.get_stencil_width() >= 1);
    // bed elevation needs more ghosts
    assert(H.get_stencil_width() >= 2);

    for (PointsWithGhosts p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // I-offset
      h_x(i, j, 0) = (H(i + 1, j) - H(i, j)) / dx;
      h_y(i, j, 0) = (+ H(i + 1, j + 1) + H(i, j + 1)
                      - H(i + 1, j - 1) - H(i, j - 1)) / (4.0*dy);
      // J-offset
      h_y(i, j, 1) = (H(i, j + 1) - H(i, j)) / dy;
      h_x(i, j, 1) = (+ H(i + 1, j + 1) + H(i + 1, j)
                      - H(i - 1, j + 1) - H(i - 1, j)) / (4.0*dx);

      // add constant bed slope in the X direction
      h_x(i, j, 0) += slope;
      h_x(i, j, 1) += slope;
      // slope in the Y direction is zero (this is here to make it
      // extra clear)
      h_y(i, j, 0) += 0.0;
      h_y(i, j, 1) += 0.0;
    }
  } else {
    SIAFD::compute_surface_gradient(h_x, h_y); 
  }

}

void SIAFD_FEvoR::compute_diffusive_flux(const IceModelVec2Stag &h_x, const IceModelVec2Stag &h_y,
                                                   IceModelVec2Stag &result, bool fast) {
  m_log->message(4,
                    "\n Starting to compute the diffusive flux\n"); 

  IceModelVec2S &thk_smooth = m_work_2d[0],
    &theta = m_work_2d[1];

  const IceModelVec2S
    &h = *m_grid->variables().get_2d_scalar("surface_altitude"),
    &H = *m_grid->variables().get_2d_scalar("land_ice_thickness");
  
  const IceModelVec2Int *mask = m_grid->variables().get_2d_mask("mask");

  bool full_update = (fast == false);

  result.set(0.0); 

  std::vector<double> delta_ij(m_grid->Mz());

  double ice_grain_size = m_config->get_double("ice_grain_size");

  bool compute_grain_size_using_age = m_config->get_boolean("compute_grain_size_using_age");

  const bool adjust_diffusivity_near_domain_boundaries = false;

  // some flow laws use grain size, and even need age to update grain size
  if (compute_grain_size_using_age && (!m_config->get_boolean("do_age"))) {
    throw RuntimeError("SIAFD::compute_diffusive_flux(): do_age not set but\n"
                       "age is needed for grain-size-based flow law");
  }

  const bool use_age = (rheology::FlowLawUsesGrainSize(flow_law()) &&
                        compute_grain_size_using_age &&
                        m_config->get_boolean("do_age"));

  // get "theta" from Schoof (2003) bed smoothness calculation and the
  // thickness relative to the smoothed bed; each IceModelVec2S involved must
  // have stencil width WIDE_GHOSTS for this too work
  m_bed_smoother->get_theta(h, theta);

  m_bed_smoother->get_smoothed_thk(h, H, *mask, thk_smooth);

  IceModelVec::AccessList list;
  list.add(theta);
  list.add(thk_smooth);
  list.add(result);

  list.add(h_x);
  list.add(h_y);

  const double *age_ij, *age_offset;
  const IceModelVec3 *age = NULL;

  if (use_age) {
    age = m_grid->variables().get_3d_scalar("age");
    list.add(*age);
  }

  if (full_update) {
    list.add(m_delta[0]);
    list.add(m_delta[1]);
  }

  const double *E_ij, *E_offset;
  const IceModelVec3 * enthalpy = m_grid->variables().get_3d_scalar("enthalpy");
  list.add(* enthalpy);


  // new code
  enhancement_factor = m_grid->variables().get_3d_scalar("enhancement_factor");
      if (enhancement_factor == NULL) {
	throw RuntimeError("enhancement_factor is not available");
      }  
  const double *EF_ij, *EF_offset;
  list.add(*enhancement_factor);
  // end of new code

  double my_D_max = 0.0;
  for (int o=0; o<2; o++) {
    for (PointsWithGhosts p(*m_grid); p; p.next()) {
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
           m_delta[o].set_column(i, j, 0.0);
        }
        continue;
      }

      if (use_age) {
	age_ij = age->get_column(i, j);
	age_offset = age->get_column(i+oi, j+oj);
      }

      E_ij = enthalpy->get_column(i, j);
      E_offset = enthalpy->get_column(i+oi, j+oj);

      // new code
      EF_ij = enhancement_factor->get_column(i, j);
      EF_offset = enhancement_factor->get_column(i+oi, j+oj);
      // end of new code

      const double slope = (o==0) ? h_x(i,j,o) : h_y(i,j,o);
      const int      ks = m_grid->kBelowHeight(thk);
      const double   alpha =
        sqrt(PetscSqr(h_x(i,j,o)) + PetscSqr(h_y(i,j,o)));
      const double theta_local = 0.5 * (theta(i,j) + theta(i+oi,j+oj));

      double  Dfoffset = 0.0;  // diffusivity for deformational SIA flow
      for (int k = 0; k <= ks; ++k) {
        double depth = thk - m_grid->z(k); // FIXME issue #15
        // pressure added by the ice (i.e. pressure difference between the
        // current level and the top of the column)
        const double pressure = m_EC->pressure(depth);

        double flow;
        if (use_age) {
          ice_grain_size = grainSizeVostok(0.5 * (age_ij[k] + age_offset[k]));
        }
        // If the flow law does not use grain size, it will just ignore it,
        // no harm there
        double E = 0.5 * (E_ij[k] + E_offset[k]);
        flow = m_flow_law->flow(alpha * pressure, E, pressure, ice_grain_size);

        // new code
        // compute the enhancement factor on the staggered grid
        double EF = 0.5 * (EF_ij[k] + EF_offset[k]);

        delta_ij[k] = EF * theta_local * 2.0 * pressure * flow;
        // end of new code

        if (k > 0) { // trapezoidal rule
          const double dz = m_grid->z(k) - m_grid->z(k-1);
          Dfoffset += 0.5 * dz * ((depth + dz) * delta_ij[k-1] + depth * delta_ij[k]);
        }
      }
      // finish off D with (1/2) dz (0 + (H-z[ks])*delta_ij[ks]), but dz=H-z[ks]:
      const double dz = thk - m_grid->z(ks);
      Dfoffset += 0.5 * dz * dz * delta_ij[ks];

      // Override diffusivity at the edges of the domain. (At these
      // locations PISM uses ghost cells *beyond* the boundary of
      // the computational domain. This does not matter if the ice
      // does not extend all the way to the domain boundary, as in
      // whole-ice-sheet simulations. In a regional setup, though,
      // this adjustment lets us avoid taking very small time-steps
      // because of the possible thickness and bed elevation
      // "discontinuities" at the boundary.)
      if (adjust_diffusivity_near_domain_boundaries && 
          (i < 0 || i >= m_grid->Mx()-1 || j < 0 || j >= m_grid->My()-1)) {
        Dfoffset = 0.0;
      }

      my_D_max = PetscMax(my_D_max, Dfoffset);

      // vertically-averaged SIA-only flux, sans sliding; note
      //   result(i,j,0) is  u  at E (east)  staggered point (i+1/2,j)
      //   result(i,j,1) is  v  at N (north) staggered point (i,j+1/2)
      result(i,j,o) = - Dfoffset * slope;

      // if doing the full update, fill the delta column above the ice and
      // store it:
      if (full_update) {
        for (unsigned int k = ks + 1; k < m_grid->Mz(); ++k) {
          delta_ij[k] = 0.0;
        }
         m_delta[o].set_column(i,j,&delta_ij[0]);
      }
    }
  } // o-loop

  m_D_max = GlobalMax(m_grid->com, my_D_max);
  m_log->message(4,
                    "\n Done computing the diffusive flux\n"); 

}

} // end of namespace stressbalance
} // end of namespace pism
