/* Copyright (C) 2023 PISM Authors
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

#ifndef PISM_GEOMETRYCALCULATOR_H
#define PISM_GEOMETRYCALCULATOR_H

#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/cell_type.hh"

namespace pism {

namespace array {
class Scalar;
}

class GeometryCalculator {
public:
  GeometryCalculator(const Config &config) {
    m_alpha = 1 - config.get_number("constants.ice.density") /
                      config.get_number("constants.sea_water.density");
    m_icefree_thickness = config.get_number("geometry.ice_free_thickness_standard");
    if (m_icefree_thickness < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid ice-free thickness threshold: %f",
                                    m_icefree_thickness);
    }
  }

  void set_icefree_thickness(double threshold) {
    if (threshold < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid ice-free thickness threshold: %f",
                                    threshold);
    }
    m_icefree_thickness = threshold;
  }

  void compute(const array::Scalar &sea_level, const array::Scalar &bed,
               const array::Scalar &thickness, array::Scalar &out_mask,
               array::Scalar &out_surface) const;

  void compute_mask(const array::Scalar &sea_level, const array::Scalar &bed,
                    const array::Scalar &thickness, array::Scalar &result) const;

  void compute_surface(const array::Scalar &sea_level, const array::Scalar &bed,
                       const array::Scalar &thickness, array::Scalar &result) const;

  inline void compute(double sea_level, double bed, double thickness, int *out_mask,
                      double *out_surface) const {
    const double hgrounded = bed + thickness; // FIXME issue #15
    const double hfloating = sea_level + m_alpha * thickness;

    const bool is_floating = (hfloating > hgrounded), ice_free = (thickness <= m_icefree_thickness);

    int mask_result;
    double surface_result;

    if (is_floating) {
      surface_result = hfloating;

      if (ice_free) {
        mask_result = cell_type::ICE_FREE_OCEAN;
      } else {
        mask_result = cell_type::ICY_OCEAN;
      }
    } else { // Grounded
      surface_result = hgrounded;

      if (ice_free) {
        mask_result = cell_type::ICE_FREE_LAND;
      } else {
        mask_result = cell_type::ICY_LAND;
      }
    }

    if (out_surface != NULL) {
      *out_surface = surface_result;
    }

    if (out_mask != NULL) {
      *out_mask = mask_result;
    }
  }

  inline int mask(double sea_level, double bed, double thickness) const {
    int result;
    compute(sea_level, bed, thickness, &result, NULL);
    return result;
  }

  inline double surface(double sea_level, double bed, double thickness) const {
    double result;
    compute(sea_level, bed, thickness, NULL, &result);
    return result;
  }

protected:
  double m_alpha;
  double m_icefree_thickness;
};

} // namespace pism

#endif /* PISM_GEOMETRYCALCULATOR_H */
