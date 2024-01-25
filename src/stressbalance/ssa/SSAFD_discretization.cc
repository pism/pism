/* Copyright (C) 2024 PISM Authors
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

#include "pism/stressbalance/ssa/SSAFD.hh"

#include "pism/util/Mask.hh"

#include "pism/rheology/FlowLaw.hh"               // rheology::averaged_hardness()
#include "pism/stressbalance/StressBalance.hh" // stressbalance::Inputs
#include "pism/geometry/Geometry.hh"           // Geometry
#include "pism/basalstrength/basal_resistance.hh" // IceBasalResistancePlasticLaw

#include "pism/util/pism_utilities.hh" // average_water_column_pressure()

namespace pism {
namespace stressbalance {

const array::Staggered & SSAFD::integrated_viscosity() const {
  return m_nuH;
}

const array::Vector& SSAFD::driving_stress() const {
  return m_taud;
}

/*!
 * Compute the weight used to determine if the difference between locations `i,j` and `n`
 * (neighbor) should be used in the computation of the surface gradient in
 * SSA::compute_driving_stress().
 *
 * We avoid differencing across
 *
 * - ice margins if stress boundary condition at ice margins (CFBC) is active
 * - grounding lines
 * - ice margins next to ice free locations above the surface elevation of the ice (fjord
 *   walls, nunataks, headwalls)
 */
static int weight(bool margin_bc,
                  int M_ij, int M_n,
                  double h_ij, double h_n,
                  int N_ij, int N_n) {
  using mask::grounded;
  using mask::icy;
  using mask::floating_ice;
  using mask::ice_free;
  using mask::ice_free_ocean;

  // grounding lines and calving fronts
  if ((grounded(M_ij) and floating_ice(M_n)) or
      (floating_ice(M_ij) and grounded(M_n)) or
      (floating_ice(M_ij) and ice_free_ocean(M_n))) {
    return 0;
  }

  // fjord walls, nunataks, headwalls
  if ((icy(M_ij) and ice_free(M_n) and h_n > h_ij) or
      (ice_free(M_ij) and icy(M_n) and h_ij > h_n)) {
    return 0;
  }

  // This condition has to match the one used to implement the calving front stress
  // boundary condition in SSAFD::assemble_rhs().
  if (margin_bc and
      ((icy(M_ij) and ice_free(M_n)) or
       (ice_free(M_ij) and icy(M_n)))) {
    return 0;
  }

  // boundaries of the "no model" area
  if ((N_ij == 0 and N_n == 1) or (N_ij == 1 and N_n == 0)) {
    return 0;
  }

  return 1;
}


//! @brief Compute the gravitational driving stress.
/*!
Computes the gravitational driving stress at the base of the ice:
\f[ \tau_d = - \rho g H \nabla h \f]
 */
void SSAFD::compute_driving_stress(const array::Scalar &ice_thickness,
                                   const array::Scalar1 &surface_elevation,
                                   const array::CellType1 &cell_type,
                                   const array::Scalar1 *no_model_mask,
                                   array::Vector &result) const {

  using mask::ice_free_ocean;
  using mask::floating_ice;

  bool cfbc = m_config->get_flag("stress_balance.calving_front_stress_bc");
  bool surface_gradient_inward =
      m_config->get_flag("stress_balance.ssa.compute_surface_gradient_inward");

  double dx = m_grid->dx(), dy = m_grid->dy();

  array::AccessScope list{ &surface_elevation, &cell_type, &ice_thickness, &result };

  if (no_model_mask != nullptr) {
    list.add(*no_model_mask);
  }

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double pressure = m_EC->pressure(ice_thickness(i, j)); // FIXME issue #15
    if (pressure <= 0.0) {
      result(i, j) = 0.0;
      continue;
    }

    // Special case for verification tests.
    if (surface_gradient_inward) {
      double h_x = diff_x_p(surface_elevation, i, j), h_y = diff_y_p(surface_elevation, i, j);
      result(i, j) = -pressure * Vector2d(h_x, h_y);
      continue;
    }

    // To compute the x-derivative we use
    //
    // * away from the grounding line, ice margins, and no_model mask transitions -- 2nd
    //   order centered difference
    //
    // * at the grounded cell near the grounding line -- 1st order
    //   one-sided difference using the grounded neighbor
    //
    // * at the floating cell near the grounding line -- 1st order
    //   one-sided difference using the floating neighbor
    //
    // All these cases can be combined by writing h_x as the weighted
    // average of one-sided differences, with weights of 0 if a finite
    // difference is not used and 1 if it is.
    //
    // The y derivative is handled the same way.

    auto M = cell_type.star_int(i, j);
    auto h = surface_elevation.star(i, j);
    stencils::Star<int> N(0);

    if (no_model_mask != nullptr) {
      N = no_model_mask->star_int(i, j);
    }

    // x-derivative
    double h_x = 0.0;
    {
      double
        west = weight(cfbc, M.c, M.w, h.c, h.w, N.c, N.w),
        east = weight(cfbc, M.c, M.e, h.c, h.e, N.c, N.e);

      if (east + west > 0) {
        h_x = 1.0 / ((west + east) * dx) * (west * (h.c - h.w) + east * (h.e - h.c));
        if (floating_ice(M.c) and (ice_free_ocean(M.e) or ice_free_ocean(M.w)))  {
          // at the ice front: use constant extrapolation to approximate the value outside
          // the ice extent (see the notes in the manual)
          h_x /= 2.0;
        }
      } else {
        h_x = 0.0;
      }
    }

    // y-derivative
    double h_y = 0.0;
    {
      double
        south = weight(cfbc, M.c, M.s, h.c, h.s, N.c, N.s),
        north = weight(cfbc, M.c, M.n, h.c, h.n, N.c, N.n);

      if (north + south > 0) {
        h_y = 1.0 / ((south + north) * dy) * (south * (h.c - h.s) + north * (h.n - h.c));
        if (floating_ice(M.c) and (ice_free_ocean(M.s) or ice_free_ocean(M.n)))  {
          // at the ice front: use constant extrapolation to approximate the value outside
          // the ice extent
          h_y /= 2.0;
        }
      } else {
        h_y = 0.0;
      }
    }

    result(i, j) = - pressure * Vector2d(h_x, h_y);
  }
}

//! \brief Checks if a cell is near or at the ice front.
/*!
 * You need to create array::AccessScope object and add `cell_type` to it.
 *
 * Note that a cell is a CFBC location of one of four direct neighbors is ice-free.
 *
 * If one of the diagonal neighbors is ice-free we don't use the CFBC, but we
 * do need to compute weights used in the SSA discretization (see
 * assemble_matrix()) to avoid differencing across interfaces between icy
 * and ice-free cells.
 *
 * This method ensures that checks in assemble_rhs() and assemble_matrix() are
 * consistent.
 */
static bool is_marginal(int i, int j, const array::CellType1 &cell_type, bool ssa_dirichlet_bc) {

  auto M = cell_type.box_int(i, j);

  using mask::ice_free;
  using mask::ice_free_ocean;
  using mask::icy;

  if (ssa_dirichlet_bc) {
    return icy(M.c) && (ice_free(M.e) || ice_free(M.w) || ice_free(M.n) || ice_free(M.s) ||
                        ice_free(M.ne) || ice_free(M.se) || ice_free(M.nw) || ice_free(M.sw));
  }

  return icy(M.c) && (ice_free_ocean(M.e) || ice_free_ocean(M.w) || ice_free_ocean(M.n) ||
                      ice_free_ocean(M.s) || ice_free_ocean(M.ne) || ice_free_ocean(M.se) ||
                      ice_free_ocean(M.nw) || ice_free_ocean(M.sw));
}

//! \brief Computes the right-hand side ("rhs") of the linear problem for the
//! Picard iteration and finite-difference implementation of the SSA equations.
/*!
The right side of the SSA equations is just the driving stress term
   \f[ - \rho g H \nabla h. \f]
The basal stress is put on the left side of the system.  This method builds the
discrete approximation of the right side.  For more about the discretization
of the SSA equations, see comments for assemble_matrix().

The values of the driving stress on the i,j grid come from a call to
compute_driving_stress().

In the case of Dirichlet boundary conditions, the entries on the right-hand side
come from known velocity values.  The fields m_bc_values and m_bc_mask are used for
this.
 */
void SSAFD::assemble_rhs(const Inputs &inputs, const array::CellType1 &cell_type,
                         const array::Vector &driving_stress,
                         array::Vector &result) {
  using mask::ice_free;
  using mask::ice_free_land;
  using mask::ice_free_ocean;

  const array::Scalar1 &bed                  = inputs.geometry->bed_elevation;
  const array::Scalar &thickness             = inputs.geometry->ice_thickness,
                      &surface               = inputs.geometry->ice_surface_elevation,
                      &sea_level             = inputs.geometry->sea_level_elevation,
                      *water_column_pressure = inputs.water_column_pressure;

  const double dx = m_grid->dx(), dy = m_grid->dy(),
               standard_gravity = m_config->get_number("constants.standard_gravity"),
               rho_ocean        = m_config->get_number("constants.sea_water.density"),
               rho_ice          = m_config->get_number("constants.ice.density");

  // This constant is for debugging: simulations should not depend on the choice of
  // velocity used in ice-free areas.
  const Vector2d ice_free_velocity(0.0, 0.0);

  const bool use_cfbc       = m_config->get_flag("stress_balance.calving_front_stress_bc"),
             flow_line_mode = m_config->get_flag("stress_balance.ssa.fd.flow_line_mode");

  // FIXME: bedrock_boundary is a misleading name
  bool bedrock_boundary = m_config->get_flag("stress_balance.ssa.dirichlet_bc");

  array::AccessScope list{ &driving_stress, &result };

  if (inputs.bc_values != nullptr and inputs.bc_mask != nullptr) {
    list.add({ inputs.bc_values, inputs.bc_mask });
  }

  if (use_cfbc) {
    list.add({ &thickness, &bed, &surface, &cell_type, &sea_level });
  }

  if (use_cfbc and (water_column_pressure != nullptr)) {
    list.add(*water_column_pressure);
  }

  result.set(0.0);

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    Vector2d taud = driving_stress(i, j);

    if (flow_line_mode) {
      // no cross-flow driving stress in the flow line mode
      taud.v = 0.0;
    }

    if ((inputs.bc_values != nullptr) and inputs.bc_mask->as_int(i, j) == 1) {
      result(i, j).u = m_scaling * (*inputs.bc_values)(i, j).u;
      result(i, j).v = m_scaling * (*inputs.bc_values)(i, j).v;
      continue;
    }

    if (use_cfbc) {
      double H_ij = thickness(i, j);

      auto M = cell_type.star_int(i, j);

      // Note: this sets velocities at both ice-free ocean and ice-free
      // bedrock to zero. This means that we need to set boundary conditions
      // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
      // to be consistent.
      if (ice_free(M.c)) {
        result(i, j) = m_scaling * ice_free_velocity;
        continue;
      }

      if (is_marginal(i, j, cell_type, bedrock_boundary)) {
        // weights at the west, east, south, and north cell faces
        int W = 0, E = 0, S = 0, N = 0;
        // direct neighbors
        // NOLINTBEGIN(readability-braces-around-statements)
        if (bedrock_boundary) {
          if (ice_free_ocean(M.e))
            E = 1;
          if (ice_free_ocean(M.w))
            W = 1;
          if (ice_free_ocean(M.n))
            N = 1;
          if (ice_free_ocean(M.s))
            S = 1;
        } else {
          if (ice_free(M.e))
            E = 1;
          if (ice_free(M.w))
            W = 1;
          if (ice_free(M.n))
            N = 1;
          if (ice_free(M.s))
            S = 1;
        }
        // NOLINTEND(readability-braces-around-statements)

        double P_ice = 0.5 * rho_ice * standard_gravity * H_ij, P_water = 0.0;

        if (water_column_pressure != nullptr) {
          P_water = (*water_column_pressure)(i, j);
        } else {
          P_water = pism::average_water_column_pressure(H_ij, bed(i, j), sea_level(i, j), rho_ice,
                                                        rho_ocean, standard_gravity);
        }

        double delta_p = H_ij * (P_ice - P_water);

        if (grid::domain_edge(*m_grid, i, j) and not(flow_line_mode or mask::grounded(M.c))) {
          // In regional setups grounded ice may extend to the edge of the domain. This
          // condition ensures that at a domain edge the ice behaves as if it extends past
          // the edge without a change in geometry.
          delta_p = 0.0;
        }

        {
          // fjord walls, nunataks, etc
          //
          // Override weights if we are at the margin and the grid cell on the other side
          // of the interface is ice-free and above the level of the ice surface.
          //
          // This effectively sets the pressure difference at the corresponding interface
          // to zero, which is exactly what we need.
          auto b   = bed.star(i, j);
          double h = surface(i, j);

          if (ice_free(M.n) and b.n > h) {
            N = 0;
          }
          if (ice_free(M.e) and b.e > h) {
            E = 0;
          }
          if (ice_free(M.s) and b.s > h) {
            S = 0;
          }
          if (ice_free(M.w) and b.w > h) {
            W = 0;
          }
        }

        // Note that if the current cell is "marginal" but not a CFBC
        // location, the following two lines are equaivalent to the "usual
        // case" below.
        //
        // Note: signs below (+E, -W, etc) are explained by directions of outward
        // normal vectors at corresponding cell faces.
        result(i, j).u = taud.u + (E - W) * delta_p / dx;
        result(i, j).v = taud.v + (N - S) * delta_p / dy;

        continue;
      } // end of "if (is_marginal(i, j))"

      // If we reached this point, then CFBC are enabled, but we are in the
      // interior of a sheet or shelf. See "usual case" below.

    } // end of "if (use_cfbc)"

    // usual case: use already computed driving stress
    result(i, j) = taud;
  }
}

//! \brief Assemble the left-hand side matrix for the KSP-based, Picard iteration,
//! and finite difference implementation of the SSA equations.
/*!
Recall the SSA equations are
\f{align*}
 - 2 \left[\nu H \left(2 u_x + v_y\right)\right]_x
        - \left[\nu H \left(u_y + v_x\right)\right]_y
        - \tau_{(b)1}  &= - \rho g H h_x, \\
   - \left[\nu H \left(u_y + v_x\right)\right]_x
        - 2 \left[\nu H \left(u_x + 2 v_y\right)\right]_y
        - \tau_{(b)2}  &= - \rho g H h_y,
\f}
where \f$u\f$ is the \f$x\f$-component of the velocity and \f$v\f$ is the
\f$y\f$-component of the velocity.

The coefficient \f$\nu\f$ is the vertically-averaged effective viscosity.
(The product \f$\nu H\f$ is computed by compute_nuH().)
The Picard iteration idea is that, to solve the nonlinear equations in which
the effective viscosity depends on the velocity, we freeze the effective
viscosity using its value at the current estimate of the velocity and we solve
the linear equations which come from this viscosity.  In abstract symbols, the
Picard iteration replaces the above nonlinear SSA equations by a sequence of
linear problems

\f[ A(U^{(k)}) U^{(k+1)} = b \f]

where \f$A(U)\f$ is the matrix from discretizing the SSA equations supposing
the viscosity is a function of known velocities \f$U\f$, and where \f$U^{(k)}\f$
denotes the \f$k\f$th iterate in the process.  The current method assembles \f$A(U)\f$.

For ice shelves \f$\tau_{(b)i} = 0\f$ [\ref MacAyealetal].
For ice streams with a basal till modelled as a plastic material,
\f$\tau_{(b)i} = - \tau_c u_i/|\mathbf{u}|\f$ where
\f$\mathbf{u} = (u,v)\f$, \f$|\mathbf{u}| = \left(u^2 + v^2\right)^{1/2}\f$,
where \f$\tau_c(t,x,y)\f$ is the yield stress of the till [\ref SchoofStream].
More generally, ice streams can be modeled with a pseudo-plastic basal till;
see IceModel::initBasalTillModel() and IceModel::updateYieldStressUsingBasalWater()
and reference [\ref BKAJS].  The pseudo-plastic till model includes all power law
sliding relations and the linearly-viscous model for sliding,
\f$\tau_{(b)i} = - \beta u_i\f$ where \f$\beta\f$ is the basal drag
(friction) parameter [\ref MacAyeal].  In any case, PISM assumes that the basal shear
stress can be factored this way, *even if the coefficient depends on the
velocity*, \f$\beta(u,v)\f$.  Such factoring is possible even in the case of
(regularized) plastic till.  This scalar coefficient \f$\beta\f$ is what is
returned by IceBasalResistancePlasticLaw::drag().

Note that the basal shear stress appears on the \em left side of the
linear system we actually solve. We believe this is crucial, because
of its effect on the spectrum of the linear approximations of each
stage. The effect on spectrum is clearest in the linearly-viscous till
case but there seems to be an analogous effect in the plastic till case.

This method assembles the matrix for the left side of the above SSA equations.
The numerical method is finite difference.  Suppose we use difference notation
\f$\delta_{+x}f^{i,j} = f^{i+1,j}-f^{i,j}\f$,
\f$\delta_{-x}f^{i,j} = f^{i,j}-f^{i-1,j}\f$, and
\f$\Delta_{x}f^{i,j} = f^{i+1,j}-f^{i-1,j}\f$, and corresponding notation for
\f$y\f$ differences, and that we write \f$N = \nu H\f$ then the first of the
two "concrete" SSA equations above has this discretization:
\f{align*}
- &2 \frac{N^{i+\frac{1}{2},j}}{\Delta x} \left[2\frac{\delta_{+x}u^{i,j}}{\Delta x} + \frac{\Delta_{y} v^{i+1,j} + \Delta_{y} v^{i,j}}{4 \Delta y}\right] + 2 \frac{N^{i-\frac{1}{2},j}}{\Delta x} \left[2\frac{\delta_{-x}u^{i,j}}{\Delta x} + \frac{\Delta_y v^{i,j} + \Delta_y v^{i-1,j}}{4 \Delta y}\right] \\
&\qquad- \frac{N^{i,j+\frac{1}{2}}}{\Delta y} \left[\frac{\delta_{+y} u^{i,j}}{\Delta y} + \frac{\Delta_x v^{i,j+1} + \Delta_x v^{i,j}}{4 \Delta x}\right] + \frac{N^{i,j-\frac{1}{2}}}{\Delta y} \left[\frac{\delta_{-y}u^{i,j}}{\Delta y} + \frac{\Delta_x v^{i,j} + \Delta_x v^{i,j-1}}{4 \Delta x}\right] - \tau_{(b)1}^{i,j} = - \rho g H^{i,j} \frac{\Delta_x h^{i,j}}{2\Delta x}.
\f}
As a picture, see Figure \ref ssastencil.

\image html ssastencil.png "\b ssastencil:  Stencil for our finite difference discretization of the first of the two scalar SSA equations.  Triangles show staggered grid points where N = nu * H is evaluated.  Circles and squares show where u and v are approximated, respectively."
\anchor ssastencil

It follows immediately that the matrix we assemble in the current method has
13 nonzeros entries per row because, for this first SSA equation, there are 5
grid values of \f$u\f$ and 8 grid values of \f$v\f$ used in this scheme.  For
the second equation we also have 13 nonzeros per row.

FIXME:  document use of DAGetMatrix and MatStencil and MatSetValuesStencil

*/
void SSAFD::fd_operator(const Inputs &inputs, const pism::Vector2d* const* input_velocity,
                        const array::Staggered1 &nuH, const array::CellType1 &cell_type, Mat *A,
                        array::Vector *Ax) {
  using mask::grounded_ice;
  using mask::ice_free;
  using mask::ice_free_land;
  using mask::ice_free_ocean;
  using mask::icy;

  const int diag_u = 4;
  const int diag_v = 13;

  PetscErrorCode ierr = 0;

  const array::Scalar1 &thickness        = inputs.geometry->ice_thickness,
                       &bed              = inputs.geometry->bed_elevation;
  const array::Scalar &surface           = inputs.geometry->ice_surface_elevation,
                      &grounded_fraction = inputs.geometry->cell_grounded_fraction,
                      &tauc              = *inputs.basal_yield_stress;

  const double dx = m_grid->dx(), dy = m_grid->dy(),
               beta_lateral_margin = m_config->get_number("basal_resistance.beta_lateral_margin"),
               beta_ice_free_bedrock =
                   m_config->get_number("basal_resistance.beta_ice_free_bedrock");

  const bool
      // FIXME: bedrock_boundary is a misleading name
      bedrock_boundary = m_config->get_flag("stress_balance.ssa.dirichlet_bc"),
      flow_line_mode   = m_config->get_flag("stress_balance.ssa.fd.flow_line_mode"),
      use_cfbc         = m_config->get_flag("stress_balance.calving_front_stress_bc"),
      replace_zero_diagonal_entries =
          m_config->get_flag("stress_balance.ssa.fd.replace_zero_diagonal_entries");

  if (A != nullptr) {
    ierr = MatZeroEntries(*A);
    PISM_CHK(ierr, "MatZeroEntries");
  }

  array::AccessScope list{ &nuH, &tauc, &cell_type, &bed, &surface };

  auto velocity = [&input_velocity](int i, int j) { return input_velocity[j][i]; };

  if (Ax != nullptr) {
    list.add(*Ax);
  }

  if (inputs.bc_values != nullptr && inputs.bc_mask != nullptr) {
    list.add(*inputs.bc_mask);
  }

  const bool sub_gl = m_config->get_flag("geometry.grounded_cell_fraction");
  if (sub_gl) {
    list.add(grounded_fraction);
  }

  // handles friction of the ice cell along ice-free bedrock margins when bedrock higher than ice
  // surface (in simplified setups)
  bool lateral_drag_enabled = m_config->get_flag("stress_balance.ssa.fd.lateral_drag.enabled");
  if (lateral_drag_enabled) {
    list.add({ &thickness, &bed, &surface });
  }
  double lateral_drag_viscosity =
      m_config->get_number("stress_balance.ssa.fd.lateral_drag.viscosity");
  double HminFrozen = 0.0;

  /* matrix assembly loop */
  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Easy cases:
      {
        // Provided Dirichlet boundary conditions
        bool bc_location = inputs.bc_mask != nullptr and inputs.bc_mask->as_int(i, j) == 1;
        // Note: this sets velocities at both ice-free ocean and ice-free
        // bedrock to zero. This means that we need to set boundary conditions
        // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
        // to be consistent.
        bool ice_free_with_cfbc = use_cfbc and cell_type.ice_free(i, j);

        if (bc_location or ice_free_with_cfbc) {
          if (Ax != nullptr) {
            (*Ax)(i, j) = m_scaling * velocity(i, j);
          }

          // set diagonal entries to one (scaled); RHS entry will be known velocity
          if (A != nullptr) {
            MatStencil row[2]  = { { -1, j, i, 0 }, { -1, j, i, 1 } };
            double identity[4] = { m_scaling, 0.0, 0.0, m_scaling };
            ierr               = MatSetValuesStencil(*A, 2, row, 2, row, identity, INSERT_VALUES);
            PISM_CHK(ierr, "MatSetValuesStencil");
          }
          continue;
        }
      }

      /* Provide shorthand for the following staggered coefficients  nu H:
       *      c_n
       *  c_w     c_e
       *      c_s
       */
      stencils::Star<double> c = nuH.star(i, j);

      if (lateral_drag_enabled) {
        // if option is set, the viscosity at ice-bedrock boundary layer will
        // be prescribed and is a temperature-independent free (user determined) parameter

        // direct neighbors
        auto M   = cell_type.star_int(i, j);
        auto H   = thickness.star(i, j);
        auto b   = bed.star(i, j);
        double h = surface(i, j);

        if (H.c > HminFrozen) {
          if (b.w > h and ice_free_land(M.w)) {
            c.w = lateral_drag_viscosity * 0.5 * (H.c + H.w);
          }
          if (b.e > h and ice_free_land(M.e)) {
            c.e = lateral_drag_viscosity * 0.5 * (H.c + H.e);
          }
          if (b.n > h and ice_free_land(M.n)) {
            c.n = lateral_drag_viscosity * 0.5 * (H.c + H.n);
          }
          if (b.s > h and ice_free_land(M.s)) {
            c.s = lateral_drag_viscosity * 0.5 * (H.c + H.s);
          }
        }
      }

      // We use DAGetMatrix to obtain the SSA matrix, which means that all 18
      // non-zeros get allocated, even though we use only 13 (or 14). The
      // remaining 5 (or 4) coefficients are zeros, but we set them anyway,
      // because this makes the code easier to understand.
      const int n_nonzeros = 18;

      // |-----+-----+---+-----+-----|
      // | NW  | NNW | N | NNE | NE  |
      // | WNW |     | | |     | ENE |
      // | W   |-----|-o-|-----| E   |
      // | WSW |     | | |     | ESE |
      // | SW  | SSW | S | SSE | SE  |
      // |-----+-----+---+-----+-----|
      //
      // We use compass rose notation for weights corresponding to interfaces between
      // cells around the current one (i, j). Here N corresponds to the interface between
      // the cell (i, j) and the one to the north of it.
      int N = 1, E = 1, S = 1, W = 1;

      // Similarly, we use compass rose notation for weights used to switch between
      // centered and one-sided finite differences. Here NNE is the interface between
      // cells N and NE, ENE - between E and NE, etc.
      int NNW = 1, NNE = 1, SSW = 1, SSE = 1;
      int WNW = 1, ENE = 1, WSW = 1, ESE = 1;

      int M_ij = cell_type.as_int(i, j);

      if (use_cfbc) {
        auto M = cell_type.box_int(i, j);

        if (is_marginal(i, j, cell_type, bedrock_boundary)) {
          // If at least one of the following four conditions is "true", we're
          // at a CFBC location.
          // NOLINTBEGIN(readability-braces-around-statements)
          if (bedrock_boundary) {

            if (ice_free_ocean(M.e))
              E = 0;
            if (ice_free_ocean(M.w))
              W = 0;
            if (ice_free_ocean(M.n))
              N = 0;
            if (ice_free_ocean(M.s))
              S = 0;

            // decide whether to use centered or one-sided differences
            if (ice_free_ocean(M.n) || ice_free_ocean(M.ne))
              NNE = 0;
            if (ice_free_ocean(M.e) || ice_free_ocean(M.ne))
              ENE = 0;
            if (ice_free_ocean(M.e) || ice_free_ocean(M.se))
              ESE = 0;
            if (ice_free_ocean(M.s) || ice_free_ocean(M.se))
              SSE = 0;
            if (ice_free_ocean(M.s) || ice_free_ocean(M.sw))
              SSW = 0;
            if (ice_free_ocean(M.w) || ice_free_ocean(M.sw))
              WSW = 0;
            if (ice_free_ocean(M.w) || ice_free_ocean(M.nw))
              WNW = 0;
            if (ice_free_ocean(M.n) || ice_free_ocean(M.nw))
              NNW = 0;

          } else { // if (not bedrock_boundary)

            if (ice_free(M.e))
              E = 0;
            if (ice_free(M.w))
              W = 0;
            if (ice_free(M.n))
              N = 0;
            if (ice_free(M.s))
              S = 0;

            // decide whether to use centered or one-sided differences
            if (ice_free(M.n) || ice_free(M.ne))
              NNE = 0;
            if (ice_free(M.e) || ice_free(M.ne))
              ENE = 0;
            if (ice_free(M.e) || ice_free(M.se))
              ESE = 0;
            if (ice_free(M.s) || ice_free(M.se))
              SSE = 0;
            if (ice_free(M.s) || ice_free(M.sw))
              SSW = 0;
            if (ice_free(M.w) || ice_free(M.sw))
              WSW = 0;
            if (ice_free(M.w) || ice_free(M.nw))
              WNW = 0;
            if (ice_free(M.n) || ice_free(M.nw))
              NNW = 0;

          } // end of the else clause following "if (bedrock_boundary)"
          // NOLINTEND(readability-braces-around-statements)
        } // end of "if (is_marginal(i, j, bedrock_boundary))"
      }   // end of "if (use_cfbc)"

      /* begin Maxima-generated code */
      const double dx2 = dx * dx, dy2 = dy * dy, d4 = 4 * dx * dy, d2 = 2 * dx * dy;

      /* Coefficients of the discretization of the first equation; u first, then v. */
      double eq1[] = {
        0,
        -c.n * N / dy2,
        0,
        -4 * c.w * W / dx2,
        (c.n * N + c.s * S) / dy2 + (4 * c.e * E + 4 * c.w * W) / dx2,
        -4 * c.e * E / dx2,
        0,
        -c.s * S / dy2,
        0,
        c.w * W * WNW / d2 + c.n * NNW * N / d4,
        (c.n * NNE * N - c.n * NNW * N) / d4 + (c.w * W * N - c.e * E * N) / d2,
        -c.e * E * ENE / d2 - c.n * NNE * N / d4,
        (c.w * W * WSW - c.w * W * WNW) / d2 + (c.n * W * N - c.s * W * S) / d4,
        (c.n * E * N - c.n * W * N - c.s * E * S + c.s * W * S) / d4 +
            (c.e * E * N - c.w * W * N - c.e * E * S + c.w * W * S) / d2,
        (c.e * E * ENE - c.e * E * ESE) / d2 + (c.s * E * S - c.n * E * N) / d4,
        -c.w * W * WSW / d2 - c.s * SSW * S / d4,
        (c.s * SSW * S - c.s * SSE * S) / d4 + (c.e * E * S - c.w * W * S) / d2,
        c.e * E * ESE / d2 + c.s * SSE * S / d4,
      };

      /* Coefficients of the discretization of the second equation; u first, then v. */
      double eq2[] = {
        c.w * W * WNW / d4 + c.n * NNW * N / d2,
        (c.n * NNE * N - c.n * NNW * N) / d2 + (c.w * W * N - c.e * E * N) / d4,
        -c.e * E * ENE / d4 - c.n * NNE * N / d2,
        (c.w * W * WSW - c.w * W * WNW) / d4 + (c.n * W * N - c.s * W * S) / d2,
        (c.n * E * N - c.n * W * N - c.s * E * S + c.s * W * S) / d2 +
            (c.e * E * N - c.w * W * N - c.e * E * S + c.w * W * S) / d4,
        (c.e * E * ENE - c.e * E * ESE) / d4 + (c.s * E * S - c.n * E * N) / d2,
        -c.w * W * WSW / d4 - c.s * SSW * S / d2,
        (c.s * SSW * S - c.s * SSE * S) / d2 + (c.e * E * S - c.w * W * S) / d4,
        c.e * E * ESE / d4 + c.s * SSE * S / d2,
        0,
        -4 * c.n * N / dy2,
        0,
        -c.w * W / dx2,
        (4 * c.n * N + 4 * c.s * S) / dy2 + (c.e * E + c.w * W) / dx2,
        -c.e * E / dx2,
        0,
        -4 * c.s * S / dy2,
        0,
      };

      /* i indices */
      const int I[] = {
        i - 1, i, i + 1, i - 1, i, i + 1, i - 1, i, i + 1,
        i - 1, i, i + 1, i - 1, i, i + 1, i - 1, i, i + 1,
      };

      /* j indices */
      const int J[] = {
        j + 1, j + 1, j + 1, j, j, j, j - 1, j - 1, j - 1,
        j + 1, j + 1, j + 1, j, j, j, j - 1, j - 1, j - 1,
      };

      /* component indices */
      const int C[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      };
      /* end Maxima-generated code */

      /* Dragging ice experiences friction at the bed determined by the
       *    IceBasalResistancePlasticLaw::drag() methods.  These may be a plastic,
       *    pseudo-plastic, or linear friction law.  Dragging is done implicitly
       *    (i.e. on left side of SSA eqns).  */
      double beta_u = 0.0, beta_v = 0.0;
      {
        double beta = 0.0;
        switch (M_ij) {
        case MASK_ICE_FREE_BEDROCK: {
          // apply drag even in this case, to help with margins; note ice free areas may
          // already have a strength extension
          beta = beta_ice_free_bedrock;
          break;
        }
        case MASK_FLOATING: {
          const Vector2d &v = velocity(i, j);
          double scaling = sub_gl ? grounded_fraction(i, j) : 0.0;
          beta = scaling * m_basal_sliding_law->drag(tauc(i, j), v.u, v.v);
          break;
        }
        case MASK_GROUNDED: {
          const Vector2d &v = velocity(i, j);
          double scaling = sub_gl ? grounded_fraction(i, j) : 1.0;
          beta = scaling * m_basal_sliding_law->drag(tauc(i, j), v.u, v.v);
          break;
        }
        case MASK_ICE_FREE_OCEAN:
        default:
          beta = 0.0;
        }

        beta_u = beta;
        beta_v = beta;
      }

      {
        // Set very high basal drag *in the direction along the boundary* at locations
        // bordering "fjord walls".

        auto M   = cell_type.star_int(i, j);
        auto b   = bed.star(i, j);
        double h = surface(i, j);

        if ((ice_free(M.n) and b.n > h) or (ice_free(M.s) and b.s > h)) {
          beta_u += beta_lateral_margin;
        }

        if ((ice_free(M.e) and b.e > h) or (ice_free(M.w) and b.w > h)) {
          beta_v += beta_lateral_margin;
        }
      }

      // add beta to diagonal entries
      eq1[diag_u] += beta_u;
      eq2[diag_v] += beta_v;

      if (flow_line_mode) {
        // set values corresponding to a trivial equation v = 0
        for (int k = 0; k < n_nonzeros; ++k) {
          eq2[k] = 0.0;
        }
        eq2[diag_v] = m_scaling;
      }

      // check diagonal entries:
      const double eps = 1e-16;
      if (fabs(eq1[diag_u]) < eps) {
        if (replace_zero_diagonal_entries) {
          eq1[diag_u] = beta_ice_free_bedrock;
        } else {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "first  (X) equation in the SSAFD system:"
                                        " zero diagonal entry at a regular (not Dirichlet B.C.)"
                                        " location: i = %d, j = %d\n",
                                        i, j);
        }
      }
      if (fabs(eq2[diag_v]) < eps) {
        if (replace_zero_diagonal_entries) {
          eq2[diag_v] = beta_ice_free_bedrock;
        } else {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "second (Y) equation in the SSAFD system:"
                                        " zero diagonal entry at a regular (not Dirichlet B.C.)"
                                        " location: i = %d, j = %d\n",
                                        i, j);
        }
      }

      // compute the matrix-vector product
      if (Ax != nullptr) {
        Vector2d sum;
        for (int k = 0; k < n_nonzeros; ++k) {
          const Vector2d &v = velocity(I[k], J[k]);
          double u_or_v = v.u * (1 - C[k]) + v.v * C[k];
          sum += Vector2d{ eq1[k], eq2[k] } * u_or_v;
        }
        (*Ax)(i, j) = sum;
      }

      // set matrix values
      if (A != nullptr) {
        MatStencil row, col[n_nonzeros];
        row.i = i;
        row.j = j;
        for (int m = 0; m < n_nonzeros; m++) {
          col[m].i = I[m];
          col[m].j = J[m];
          col[m].c = C[m];
        }

        // set coefficients of the first equation:
        row.c = 0;
        ierr  = MatSetValuesStencil(*A, 1, &row, n_nonzeros, col, eq1, INSERT_VALUES);
        PISM_CHK(ierr, "MatSetValuesStencil");

        // set coefficients of the second equation:
        row.c = 1;
        ierr  = MatSetValuesStencil(*A, 1, &row, n_nonzeros, col, eq2, INSERT_VALUES);
        PISM_CHK(ierr, "MatSetValuesStencil");
      }
    } // i,j-loop
  } catch (...) {
    loop.failed();
  }
  loop.check();

  if (A != nullptr) {
    ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
    PISM_CHK(ierr, "MatAssemblyBegin");

    ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
    PISM_CHK(ierr, "MatAssemblyEnd");
#if (Pism_DEBUG == 1)
    ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    PISM_CHK(ierr, "MatSetOption");
#endif
  }
}

//! @brief Computes vertically-averaged ice hardness on the staggered grid.
void SSAFD::compute_average_ice_hardness(const array::Scalar1 &ice_thickness,
                                         const array::Array3D &enthalpy,
                                         const array::CellType1 &cell_type,
                                         array::Staggered &result) {

  const double *E_ij = nullptr, *E_offset = nullptr;

  auto Mz = m_grid->Mz();
  std::vector<double> E(Mz);

  array::AccessScope list{ &ice_thickness, &enthalpy, &result, &cell_type };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      E_ij = enthalpy.get_column(i, j);
      for (int o = 0; o < 2; o++) {
        const int oi = 1 - o, oj = o;
        double H;

        if (cell_type.icy(i, j) && cell_type.icy(i + oi, j + oj)) {
          H = 0.5 * (ice_thickness(i, j) + ice_thickness(i + oi, j + oj));
        } else if (cell_type.icy(i, j)) {
          H = ice_thickness(i, j);
        } else {
          H = ice_thickness(i + oi, j + oj);
        }

        if (H == 0) {
          result(i, j, o) = -1e6; // an obviously impossible value
          continue;
        }

        E_offset = enthalpy.get_column(i + oi, j + oj);
        // build a column of enthalpy values a the current location:
        for (unsigned int k = 0; k < Mz; ++k) {
          E[k] = 0.5 * (E_ij[k] + E_offset[k]);
        }

        result(i, j, o) = rheology::averaged_hardness(*m_flow_law, H, m_grid->kBelowHeight(H),
                                                      m_grid->z().data(), E.data());
      } // o
    }   // loop over points
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

/*!
 *  These computations do not depend on the solution, so they need to be done only once.
 *
 *  Updates m_cell_type, m_taud, m_hardness.
 */
void SSAFD::initialize_iterations(const Inputs &inputs) {
  // update the cell type mask using the ice-free thickness threshold for stress balance
  // computations
  {
    const double H_threshold = m_config->get_number("stress_balance.ice_free_thickness_standard");
    GeometryCalculator gc(*m_config);
    gc.set_icefree_thickness(H_threshold);

    gc.compute_mask(inputs.geometry->sea_level_elevation, inputs.geometry->bed_elevation,
                    inputs.geometry->ice_thickness, //
                    m_cell_type);                   // output
    // note: compute_mask() updates ghosts without communication ("redundantly")
  }
  compute_driving_stress(inputs.geometry->ice_thickness, inputs.geometry->ice_surface_elevation,
                         m_cell_type, inputs.no_model_mask, //
                         m_taud);                           // output
  compute_average_ice_hardness(inputs.geometry->ice_thickness,
                               *inputs.enthalpy,
                               m_cell_type, m_hardness);

  if (inputs.fracture_density != nullptr) {
    fracture_induced_softening(*inputs.fracture_density, m_hardness);
  }
}

/*! @brief Correct vertically-averaged hardness using a
    parameterization of the fracture-induced softening.

  See T. Albrecht, A. Levermann; Fracture-induced softening for
  large-scale ice dynamics; (2013), The Cryosphere Discussions 7;
  4501-4544; DOI:10.5194/tcd-7-4501-2013

  Note that this paper proposes an adjustment of the enhancement factor:

  \f[E_{\text{effective}} = E \cdot (1 - (1-\epsilon) \phi)^{-n}.\f]

  Let \f$E_{\text{effective}} = E\cdot C\f$, where \f$C\f$ is the
  factor defined by the formula above.

  Recall that the effective viscosity is defined by

  \f[\nu(D) = \frac12 B D^{(1-n)/(2n)}\f]

  and the viscosity form of the flow law is

  \f[\sigma'_{ij} = E_{\text{effective}}^{-\frac1n}2\nu(D) D_{ij}.\f]

  Then

  \f[\sigma'_{ij} = E_{\text{effective}}^{-\frac1n}BD^{(1-n)/(2n)}D_{ij}.\f]

  Using the fact that \f$E_{\text{effective}} = E\cdot C\f$, this can be rewritten as

  \f[\sigma'_{ij} = E^{-\frac1n} \left(C^{-\frac1n}B\right) D^{(1-n)/(2n)}D_{ij}.\f]

  So scaling the enhancement factor by \f$C\f$ is equivalent to scaling
  ice hardness \f$B\f$ by \f$C^{-\frac1n}\f$.
*/
void SSAFD::fracture_induced_softening(const array::Scalar &fracture_density,
                                       array::Staggered &ice_hardness) {

  const double epsilon = m_config->get_number("fracture_density.softening_lower_limit"),
               n_glen  = m_flow_law->exponent();

  array::AccessScope list{ &ice_hardness, &fracture_density };

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (int o = 0; o < 2; o++) {
      const int oi = 1 - o, oj = o;

      const double
          // fracture density on the staggered grid:
          phi = 0.5 * (fracture_density(i, j) + fracture_density(i + oi, j + oj)),
          // the line below implements equation (6) in the paper
          softening = pow((1.0 - (1.0 - epsilon) * phi), -n_glen);

      ice_hardness(i, j, o) *= pow(softening, -1.0 / n_glen);
    }
  }
}

//! \brief Compute the product of ice thickness and effective viscosity (on the
//! staggered grid).
/*!
In PISM the product \f$\nu H\f$ can be
  - constant, or
  - can be computed with a constant ice hardness \f$\bar B\f$ (temperature-independent)
    but with dependence of the viscosity on the strain rates, or
  - it can depend on the strain rates \e and have a vertically-averaged ice
    hardness depending on temperature or enthalpy.

The flow law in ice stream and ice shelf regions must, for now, be a
(temperature-dependent) Glen law.  This is the only flow law we know how to
invert to ``viscosity form''.  More general forms like Goldsby-Kohlstedt are
not yet inverted.

The viscosity form of a Glen law is
   \f[ \nu(T^*,D) = \frac{1}{2} B(T^*) D^{(1/n)-1}\, D_{ij} \f]
where
   \f[  D_{ij} = \frac{1}{2} \left(\frac{\partial U_i}{\partial x_j} +
                                   \frac{\partial U_j}{\partial x_i}\right) \f]
is the strain rate tensor and \f$B\f$ is an ice hardness related to
the ice softness \f$A(T^*)\f$ by
   \f[ B(T^*)=A(T^*)^{-1/n}  \f]
in the case of a temperature dependent Glen-type law.  (Here \f$T^*\f$ is the
pressure-adjusted temperature.)

The effective viscosity is then
   \f[ \nu = \frac{\bar B}{2} \left[\left(\frac{\partial u}{\partial x}\right)^2 +
                               \left(\frac{\partial v}{\partial y}\right)^2 +
                               \frac{\partial u}{\partial x} \frac{\partial v}{\partial y} +
                               \frac{1}{4} \left(\frac{\partial u}{\partial y}
                                                 + \frac{\partial v}{\partial x}\right)^2
                               \right]^{(1-n)/(2n)}  \f]
where in the temperature-dependent case
   \f[ \bar B = \frac{1}{H}\,\int_b^h B(T^*)\,dz\f]
This integral is approximately computed by the trapezoid rule.

The result of this procedure is \f$\nu H\f$, not just \f$\nu\f$, this it is
a vertical integral, not average, of viscosity.

The resulting effective viscosity times thickness is regularized by ensuring that
its minimum is at least \f$\epsilon\f$.  This regularization constant is an argument.

In this implementation we set \f$\nu H\f$ to a constant anywhere the ice is
thinner than a certain minimum. See SSAStrengthExtension and compare how this
issue is handled when -cfbc is set.
*/
void SSAFD::compute_nuH_everywhere(const array::Scalar1 &ice_thickness,
                                   const pism::Vector2d *const *velocity,
                                   const array::Staggered &hardness, double nuH_regularization,
                                   array::Staggered &result) {

  auto uv = [&velocity](int i, int j) { return velocity[j][i]; };

  array::AccessScope list{ &result, &hardness, &ice_thickness };

  double n_glen                 = m_flow_law->exponent(),
         nu_enhancement_scaling = 1.0 / pow(m_e_factor, 1.0 / n_glen);

  const double dx = m_grid->dx(), dy = m_grid->dy();

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    const int
      E = i + 1,
      W = i - 1,
      N = j + 1,
      S = j - 1;

    stencils::Box<Vector2d> V{ uv(i, j), uv(i, N), uv(W, N), uv(W, j), uv(W, S),
                               uv(i, S), uv(E, S), uv(E, j), uv(E, N) };

    for (int o = 0; o < 2; ++o) {
      const int oi = 1 - o, oj = o;

      const double H = 0.5 * (ice_thickness(i, j) + ice_thickness(i + oi, j + oj));

      if (H < strength_extension->get_min_thickness()) {
        result(i, j, o) = strength_extension->get_notional_strength();
        continue;
      }

      double u_x, u_y, v_x, v_y;
      // Check the offset to determine how to differentiate velocity
      if (o == 0) {
        u_x = (V.e.u - V.c.u) / dx;
        v_x = (V.e.v - V.c.v) / dx;
        u_y = (V.n.u + V.ne.u - V.s.u - V.se.u) / (4 * dy);
        v_y = (V.n.v + V.ne.v - V.s.v - V.se.v) / (4 * dy);
      } else {
        u_x = (V.e.u + V.ne.u - V.w.u - V.nw.u) / (4 * dx);
        v_x = (V.e.v + V.ne.v - V.w.v - V.nw.v) / (4 * dx);
        u_y = (V.n.u - V.c.u) / dy;
        v_y = (V.n.v - V.c.v) / dy;
      }

      double nu = 0.0;
      m_flow_law->effective_viscosity(hardness(i, j, o),
                                      secondInvariant_2D({ u_x, v_x }, { u_y, v_y }), &nu, nullptr);

      result(i, j, o) = nu * H;

      // include the SSA enhancement factor; in most cases m_e_factor is 1
      result(i, j, o) *= nu_enhancement_scaling;

      // We ensure that nuH is bounded below by a positive constant.
      result(i, j, o) += nuH_regularization;
    } // o-loop
  }   // i,j-loop
}

/**
 * @brief Compute the product of ice viscosity and thickness on the
 * staggered grid. Used when CFBC is enabled.
 *
 * 1) Loops over all grid points and width=1 ghosts and estimates u_x and v_x at the
 *    i-offset staggered grid locations and u_y and v_y at the j-offset staggered grid
 *    locations. This requires width=2 ghosts of `velocity` and `cell_type`.
 *
 * 2) Loops over all grid points (excluding ghost points) and computes weighted averages
 *    of quantities from step 1 to estimate (u_y, v_y) at i-offset locations and (u_x,
 *    v_x) and j-offset locations. This uses width=1 ghost values set in step 1.
 *
 * 3) In the second loop, ice thickness, ice hardness and (u_x, u_y, v_x, v_y) at
 *    staggered grid locations from steps 1 and 2 are used to estimate nuH.
 *
 * @param[out] result nu*H product
 * @param[in] nuH_regularization regularization parameter (added to nu*H to keep it away from zero)
 *
 * @return 0 on success
 */
void SSAFD::compute_nuH_cfbc(const array::Scalar1 &ice_thickness, const array::CellType2 &cell_type,
                             const pism::Vector2d* const* velocity,
                             const array::Staggered &hardness,
                             double nuH_regularization, array::Staggered &result) {

  double n_glen                 = m_flow_law->exponent(),
         nu_enhancement_scaling = 1.0 / pow(m_e_factor, 1.0 / n_glen);

  const double dx = m_grid->dx(), dy = m_grid->dy();

  array::AccessScope list{ &cell_type, &m_work };

  assert(m_work.stencil_width() >= 1);

  auto U = [&velocity](int i, int j) { return velocity[j][i].u; };
  auto V = [&velocity](int i, int j) { return velocity[j][i].v; };

  for (auto p = m_grid->points(1); p; p.next()) {
    const int i = p.i(), j = p.j();

    // x-derivative, i-offset
    {
      if (cell_type.icy(i, j) && cell_type.icy(i + 1, j)) {
        m_work(i, j).u_x = (U(i + 1, j) - U(i, j)) / dx; // u_x
        m_work(i, j).v_x = (V(i + 1, j) - V(i, j)) / dx; // v_x
        m_work(i, j).w_i = 1.0;
      } else {
        m_work(i, j).u_x = 0.0;
        m_work(i, j).v_x = 0.0;
        m_work(i, j).w_i = 0.0;
      }
    }

    // y-derivative, j-offset
    {
      if (cell_type.icy(i, j) && cell_type.icy(i, j + 1)) {
        m_work(i, j).u_y = (U(i, j + 1) - U(i, j)) / dy; // u_y
        m_work(i, j).v_y = (V(i, j + 1) - V(i, j)) / dy; // v_y
        m_work(i, j).w_j = 1.0;
      } else {
        m_work(i, j).u_y = 0.0;
        m_work(i, j).v_y = 0.0;
        m_work(i, j).w_j = 0.0;
      }
    }
  } // end of the loop over grid points and width=1 ghosts

  list.add({ &result, &hardness, &ice_thickness });

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double u_x, u_y, v_x, v_y, H, nu, W;
    // i-offset
    {
      if (cell_type.icy(i, j) && cell_type.icy(i + 1, j)) {
        H = 0.5 * (ice_thickness(i, j) + ice_thickness(i + 1, j));
      } else if (cell_type.icy(i, j)) {
        H = ice_thickness(i, j);
      } else {
        H = ice_thickness(i + 1, j);
      }

      if (H >= strength_extension->get_min_thickness()) {
        u_x = m_work(i, j).u_x;
        v_x = m_work(i, j).v_x;

        W = m_work(i, j).w_j + m_work(i, j - 1).w_j + m_work(i + 1, j - 1).w_j +
            m_work(i + 1, j).w_j;
        if (W > 0) {
          u_y = 1.0 / W *
                (m_work(i, j).u_y + m_work(i, j - 1).u_y + m_work(i + 1, j - 1).u_y +
                 m_work(i + 1, j).u_y);
          v_y = 1.0 / W *
                (m_work(i, j).v_y + m_work(i, j - 1).v_y + m_work(i + 1, j - 1).v_y +
                 m_work(i + 1, j).v_y);
        } else {
          u_y = 0.0;
          v_y = 0.0;
        }

        m_flow_law->effective_viscosity(
            hardness(i, j, 0), secondInvariant_2D({ u_x, v_x }, { u_y, v_y }), &nu, nullptr);
        result(i, j, 0) = nu * H;
      } else {
        result(i, j, 0) = strength_extension->get_notional_strength();
      }
    }

    // j-offset
    {
      if (cell_type.icy(i, j) && cell_type.icy(i, j + 1)) {
        H = 0.5 * (ice_thickness(i, j) + ice_thickness(i, j + 1));
      } else if (cell_type.icy(i, j)) {
        H = ice_thickness(i, j);
      } else {
        H = ice_thickness(i, j + 1);
      }

      if (H >= strength_extension->get_min_thickness()) {
        u_y = m_work(i, j).u_y;
        v_y = m_work(i, j).v_y;

        W = m_work(i, j).w_i + m_work(i - 1, j).w_i + m_work(i - 1, j + 1).w_i +
            m_work(i, j + 1).w_i;
        if (W > 0.0) {
          u_x = 1.0 / W *
                (m_work(i, j).u_x + m_work(i - 1, j).u_x + m_work(i - 1, j + 1).u_x +
                 m_work(i, j + 1).u_x);
          v_x = 1.0 / W *
                (m_work(i, j).v_x + m_work(i - 1, j).v_x + m_work(i - 1, j + 1).v_x +
                 m_work(i, j + 1).v_x);
        } else {
          u_x = 0.0;
          v_x = 0.0;
        }

        m_flow_law->effective_viscosity(hardness(i, j, 1),
                                        secondInvariant_2D({ u_x, v_x }, { u_y, v_y }), &nu, NULL);
        result(i, j, 1) = nu * H;
      } else {
        result(i, j, 1) = strength_extension->get_notional_strength();
      }
    }

    // adjustments:
    for (int o = 0; o < 2; ++o) {
      // include the SSA enhancement factor; in most cases ssa_enhancement_factor is 1
      result(i, j, o) *= nu_enhancement_scaling;

      // We ensure that nuH is bounded below by a positive constant.
      result(i, j, o) += nuH_regularization;
    }
  } // end of the loop over grid points
}

void SSAFD::compute_nuH(const array::Scalar1 &ice_thickness, const array::CellType2 &cell_type,
                        const pism::Vector2d *const *velocity, const array::Staggered &hardness,
                        double nuH_regularization, array::Staggered1 &result) {
  if (m_config->get_flag("stress_balance.calving_front_stress_bc")) {
    compute_nuH_cfbc(ice_thickness, cell_type, velocity, hardness, nuH_regularization, result);
  } else {
    compute_nuH_everywhere(ice_thickness, velocity, m_hardness, nuH_regularization, result);
  }
  result.update_ghosts();
}

void SSAFD::compute_residual(const Inputs &inputs, const array::Vector &velocity,
                             array::Vector &result) {
  // update m_cell_type, m_taud, m_hardness
  initialize_iterations(inputs);

  m_velocity.copy_from(velocity);

  {
    array::AccessScope list{ &m_velocity };
    compute_nuH(inputs.geometry->ice_thickness, m_cell_type, m_velocity.array(), m_hardness,
                m_config->get_number("stress_balance.ssa.epsilon"), m_nuH);

    fd_operator(inputs, m_velocity.array(), m_nuH, m_cell_type, nullptr, &result);
  }
  assemble_rhs(inputs, m_cell_type, m_taud, m_rhs);

  result.add(-1.0, m_rhs);
}

} // namespace stressbalance
} // namespace pism
