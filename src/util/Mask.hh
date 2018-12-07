// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Constantine Khroulev and David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef _MASK_H_
#define _MASK_H_

// the following three includes are needed here because of inlined code
#include "iceModelVec.hh"
#include "ConfigInterface.hh"
#include "error_handling.hh"

namespace pism {

enum MaskValue {
  MASK_UNKNOWN          = -1,
  MASK_ICE_FREE_BEDROCK = 0,
  MASK_GROUNDED         = 2,
  MASK_FLOATING         = 3,
  MASK_ICE_FREE_OCEAN   = 4
};

namespace mask {
//! \brief An ocean cell (floating ice or ice-free).
  inline bool ocean(int M) {
    return M >= MASK_FLOATING;
  }
  //! \brief Grounded cell (grounded ice or ice-free).
  inline bool grounded(int M) {
    return not ocean(M);
  }
  //! \brief Ice-filled cell (grounded or floating).
  inline bool icy(int M) {
    return (M == MASK_GROUNDED) || (M == MASK_FLOATING);
  }
  inline bool grounded_ice(int M) {
    return icy(M) && grounded(M);
  }
  inline bool floating_ice(int M) {
    return icy(M) && ocean(M);
  }
  //! \brief Ice-free cell (grounded or ocean).
  inline bool ice_free(int M) {
    return not icy(M);
  }
  inline bool ice_free_ocean(int M) {
    return ocean(M) && ice_free(M);
  }
  inline bool ice_free_land(int M) {
    return grounded(M) && ice_free(M);
  }
}

class GeometryCalculator {
public:
  GeometryCalculator(const Config &config) {
    m_alpha = 1 - config.get_double("constants.ice.density") / config.get_double("constants.sea_water.density");
    m_is_dry_simulation = config.get_boolean("ocean.always_grounded");
    m_icefree_thickness = config.get_double("geometry.ice_free_thickness_standard");
    if (m_icefree_thickness < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "invalid ice-free thickness threshold: %f", m_icefree_thickness);
    }
  }

  void set_icefree_thickness(double threshold) {
    if (threshold < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid ice-free thickness threshold: %f", threshold);
    }
    m_icefree_thickness = threshold;
  }

  void compute(const IceModelVec2S &sea_level, const IceModelVec2S &bed, const IceModelVec2S &thickness,
               IceModelVec2Int &out_mask, IceModelVec2S &out_surface) const;

  void compute_mask(const IceModelVec2S& sea_level, const IceModelVec2S& bed,
                    const IceModelVec2S& thickness, IceModelVec2Int& result) const;

  void compute_surface(const IceModelVec2S& sea_level, const IceModelVec2S& bed,
                       const IceModelVec2S& thickness, IceModelVec2S& result) const;

  inline void compute(double sea_level, double bed, double thickness,
                      int *out_mask, double *out_surface) const {
    const double hgrounded = bed + thickness; // FIXME issue #15
    const double hfloating = sea_level + m_alpha*thickness;

    const bool
      is_floating = (hfloating > hgrounded),
      ice_free    = (thickness <= m_icefree_thickness);

    int mask_result;
    double surface_result;

    if (is_floating && (not m_is_dry_simulation)) {
      surface_result = hfloating;

      if (ice_free) {
        mask_result = MASK_ICE_FREE_OCEAN;
      } else {
        mask_result = MASK_FLOATING;
      }
    } else {  // Grounded
      surface_result = hgrounded;

      if (ice_free) {
        mask_result = MASK_ICE_FREE_BEDROCK;
      } else {
        mask_result = MASK_GROUNDED;
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
  bool m_is_dry_simulation;
};

} // end of namespace pism

#endif /* _MASK_H_ */
