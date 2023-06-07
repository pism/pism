/* Copyright (C) 2016, 2019, 2020, 2022, 2023 PISM Authors
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

#ifndef PISM_ARRAY_CELLTYPE_H
#define PISM_ARRAY_CELLTYPE_H

#include "pism/util/array/Scalar.hh"
#include "pism/util/cell_type.hh"

namespace pism {
namespace array {

//! "Cell type" mask. Adds convenience methods to `array::Scalar`.
class CellType : public Scalar {
public:
  CellType(std::shared_ptr<const Grid> grid, const std::string &name);

  inline bool water(int i, int j) const {
    return cell_type::water(as_int(i, j));
  }

  inline bool land(int i, int j) const {
    return cell_type::land(as_int(i, j));
  }

  inline bool icy(int i, int j) const {
    return cell_type::icy(as_int(i, j));
  }

  inline bool grounded_ice(int i, int j) const {
    return cell_type::grounded_ice(as_int(i, j));
  }

  inline bool floating_ice(int i, int j) const {
    return cell_type::floating_ice(as_int(i, j));
  }

  inline bool ice_free(int i, int j) const {
    return cell_type::ice_free(as_int(i, j));
  }

  // inline bool ice_free_ocean(int i, int j) const {
  //   return cell_type::ice_free_ocean(as_int(i, j));
  // }

  inline bool ice_free_water(int i, int j) const {
    return cell_type::ice_free_water(as_int(i, j));
  }

  inline bool ice_free_land(int i, int j) const {
    return cell_type::ice_free_land(as_int(i, j));
  }
protected:
  CellType(std::shared_ptr<const Grid> grid, const std::string &name, int w);
};

/*!
 * Cell type array supporting width=1 stencil computations (ghosted).
 */
class CellType1 : public CellType {
public:
  CellType1(std::shared_ptr<const Grid> grid, const std::string &name);
  using Array2D<double>::star;
  using Array2D<double>::box;
  using Scalar::star_int;
  using Scalar::box_int;

  //! \brief Ice margin (ice-filled with at least one of four neighbors ice-free).
  inline bool ice_margin(int i, int j) const {
    return icy(i, j) and (ice_free(i + 1, j) or ice_free(i - 1, j) or
                          ice_free(i, j + 1) or ice_free(i, j - 1));
  }

  //! \brief Ice-free margin (at least one of four neighbors has ice).
  inline bool next_to_ice(int i, int j) const {
    return (icy(i + 1, j) or icy(i - 1, j) or icy(i, j + 1) or icy(i, j - 1));
  }

  inline bool next_to_floating_ice(int i, int j) const {
    return (floating_ice(i + 1, j) or floating_ice(i - 1, j) or
            floating_ice(i, j + 1) or floating_ice(i, j - 1));
  }

  inline bool next_to_grounded_ice(int i, int j) const {
    return (grounded_ice(i + 1, j) or grounded_ice(i - 1, j) or
            grounded_ice(i, j + 1) or grounded_ice(i, j - 1));
  }

  inline bool next_to_ice_free_land(int i, int j) const {
    return (ice_free_land(i + 1, j) or ice_free_land(i - 1, j) or
            ice_free_land(i, j + 1) or ice_free_land(i, j - 1));
  }

  // inline bool next_to_ice_free_ocean(int i, int j) const {
  //   return (ice_free_ocean(i + 1, j) or ice_free_ocean(i - 1, j) or
  //           ice_free_ocean(i, j + 1) or ice_free_ocean(i, j - 1));
  // }

  inline bool next_to_ice_free_water(int i, int j) const {
    return (ice_free_water(i + 1, j) or ice_free_water(i - 1, j) or
            ice_free_water(i, j + 1) or ice_free_water(i, j - 1));
  }
protected:
  CellType1(std::shared_ptr<const Grid> grid, const std::string &name, int width);
};

/*!
 * Cell type array supporting width=2 stencil computations (ghosted).
 */
class CellType2 : public CellType1 {
public:
  CellType2(std::shared_ptr<const Grid> grid, const std::string &name);
};

} // end of namespace array
} // end of namespace pism

#endif /* PISM_ARRAY_CELLTYPE_H */
