/* Copyright (C) 2016 PISM Authors
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

#include "base/util/iceModelVec.hh"
#include "base/util/error_handling.hh"
#include "base/util/Mask.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {

/**
 * Compute `lambda`, the grounding line position (LI in [@ref Gladstoneetal2010]).
 *
 * @param mu `mu = ice_density / water_density`
 * @param sea_level sea level elevation
 * @param H0 ice thickness at the first point
 * @param b0 bed elevation at the first point
 * @param H1 ice thickness at the second point
 * @param b1 bed elevation at the second point
 *
 * Here `mu * H` is the depth of the ice below sea water and `sea_level - b` is the ocean depth.
 *
 * If `mu * H` is greater than `sea_level - b` then there is too much ice for ice to float.
 *
 * @return `lambda`, a number between 0 and 1
 */
static inline double gl_position(double mu,
                                 double sea_level,
                                 double H0, double b0,
                                 double H1, double b1) {
  const double
    alpha = mu * H0 - (sea_level - b0),
    beta = mu * H1 - (sea_level - b1);

  if (alpha - beta == 0.0) {
    throw RuntimeError("cannot determine grounding line position. Please submit a bug report.");
  }

  double lambda = alpha / (alpha - beta);

  return std::max(0.0, std::min(lambda, 1.0));
}


/**
   @brief Updates the fractional "flotation mask".

   This mask ranges from 0 to 1 and is equal to the fraction of the
   cell (by area) that is grounded.

   Currently it is used to adjust the basal drag near the grounding
   line in the SSAFD stress balance model and the basal melt rate
   computation in the energy balance code (IceModel::temperatureStep()
   and IceModel::enthalpyAndDrainageStep()).

   We use the 1D (flow line) parameterization of the sub-grid
   grounding line position due to [@ref Gladstoneetal2010], (section
   3.1.1) and generalize it to the case of arbitrary sea level
   elevation. Then this sub-grid grounding line position is used to
   compute the grounded area fraction for each cell.

   We use the "LI" 1D grounding line position parameterization in "x"
   and "y" directions *independently*. The grounding line positions in
   the "x" and "y" directions (@f$ \lambda_x @f$ and @f$ \lambda_y
   @f$) are then interpreted as *width* and *height* of a rectangular
   sub-set of the cell. The grounded area is computed as the product
   of these two: @f$ A_g = \lambda_x \times \lambda_y. @f$

   Consider a cell at `(i,j)` and assume that the ice is grounded
   there and floating at `(i+1,j)`.

   Assume that the ice thickness and bedrock elevation change linearly
   from `(i,j)` to `(i+1,j)` and drop the `j` index for clarity:

   @f{align*}{
   H(\lambda) &= H_{i}(1 - \lambda) + H_{i+1}\lambda,\\
   b(\lambda) &= b_{i}(1 - \lambda) + b_{i+1}\lambda.\\
   @f}

   Here @f$ \lambda @f$ is the dimensionless variable parameterizing
   the sub-grid grounding line position, ranging from 0 to 1.

   Now, substituting @f$ b(\lambda) @f$ and @f$ H(\lambda) @f$ into
   the floatation criterion

   @f{align*}{
   \mu\cdot H(\lambda) &= z_{\text{sea level}} - b(\lambda),\\
   \mu &= \rho_{\text{ice}} / \rho_{\text{sea water}}.
   @f}

   and solving for @f$ \lambda @f$, we get

   @f{align*}{
   \lambda_{g} &= \frac{\alpha}{\alpha - \beta}\\
   & \text{where} \\
   \alpha &= \mu\cdot H_{i} + b_{i} - z_{\text{sea level}}\\
   \beta &= \mu\cdot H_{i+1} + b_{i+1} - z_{\text{sea level}}.
   @f}

   Note that [@ref Gladstoneetal2010] describe a parameterization of
   the grounding line position within a cell defined as the interval
   from the grid point `(i)` to the grid point `(i+1)`, with the ice
   thickness @f$ H @f$ and the bed elevation @f$ b @f$ defined *at
   grid points* (boundaries of a 1D cell).

   Here we compute a grounded fraction of the **cell centered at the
   grid point**.

   FIXME: sometimes alpha<0 (slightly below flotation) even though the
   mask says grounded

   @param[in] ice_density ice density, kg/m3
   @param[in] ocean_density ocean (sea water) density, kg/m3
   @param[in] sea_level sea level elevation, m
   @param[in] ice_thickness ice thickness (a 2D field)
   @param[in] bed_topography bed topography (a 2D field)
   @param[in] mask cell type mask (a 2D field)
   @param[out] result resulting grounded fraction
   @param[out] result_x grounding line position in the x direction (1D parameterization, for debugging)
   @param[out] result_y grounding line position in the y direction (1D parameterization, for debugging)
 */
void compute_grounded_cell_fraction(double ice_density,
                             double ocean_density,
                             double sea_level,
                             const IceModelVec2S &ice_thickness,
                             const IceModelVec2S &bed_topography,
                             const IceModelVec2CellType &mask,
                             IceModelVec2S &result,
                             IceModelVec2S *result_x,
                             IceModelVec2S *result_y) {

  const double mu = ice_density / ocean_density;

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(bed_topography);
  list.add(mask);
  list.add(result);

  if (result_x != NULL) {
    list.add(*result_x);
  }
  if (result_y != NULL) {
    list.add(*result_y);
  }

  IceGrid::ConstPtr grid = result.get_grid();

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const StarStencil<int>
        m = mask.int_star(i, j);
      const StarStencil<double>
        H = ice_thickness.star(i, j),
        b = bed_topography.star(i, j);

      const Direction dirs[4] = {North, East, South, West};

      // Contributions from the 4 directions
      std::vector<double> lambda(4, 0.0);

      if (mask::grounded(m.ij)) {

        for (int k = 0; k < 4; ++k) {
          const Direction direction = dirs[k];

          if (mask::grounded(m[direction])) {
            lambda[k] = 0.5;
          } else {
            const double L = gl_position(mu, sea_level, H.ij, b.ij, H[direction], b[direction]);
            lambda[k] = std::min(L, 0.5);
          }
        }

      } else if (mask::ocean(m.ij)) {

        for (int k = 0; k < 4; ++k) {
          const Direction direction = dirs[k];

          if (mask::grounded(m[direction])) {
            const double L = gl_position(mu, sea_level, H[direction], b[direction], H.ij, b.ij);
            lambda[k] = std::max(L - 0.5, 0.0);
          } else {
            lambda[k] = 0.0;
          }
        }

      } else { // end of the "mask::ocean(m.ij)" case
        throw RuntimeError::formatted("ice is neither grounded nor floating (ocean) at (%d,%d)."
                                      " This should not happen.", i, j);
      }

      const double
        lambda_x = lambda[East] + lambda[West],
        lambda_y = lambda[South] + lambda[North];

      if (result_x != NULL) {
        (*result_x)(i, j) = lambda_x;
      }
      if (result_y != NULL) {
        (*result_y)(i, j) = lambda_y;
      }

      // Note: in the flowline case (in general: when one of lambda_[xy] is zero) this does not work
      // very well, but I don't know how to get it "right"... (CK)
      result(i, j) = lambda_x * lambda_y;

    } // end of the loop over grid points
  } catch (...) {
    loop.failed();
  }
  loop.check();

}

} // end of namespace pism
