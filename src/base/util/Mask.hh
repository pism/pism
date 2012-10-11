// Copyright (C) 2011, 2012 Constantine Khroulev and David Maxwell
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

#ifndef _MASK_H_
#define _MASK_H_

// the following three includes are needed here because of inlined code
#include "iceModelVec.hh"
#include "NCVariable.hh"
#include "flowlaws.hh"
#include "IceGrid.hh"

class Mask
{
public:
  Mask() {}
  ~Mask() {}
  //! \brief An ocean cell (floating ice or ice-free).
  inline bool ocean(int M) { return M >= MASK_FLOATING; }
  //! \brief Grounded cell (grounded ice or ice-free).
  inline bool grounded(int M) { return !ocean(M); }
  //! \brief Ice-filled cell (grounded or floating).
  inline bool icy(int M) { return (M == MASK_GROUNDED) || (M == MASK_FLOATING); }
  inline bool grounded_ice(int M) { return icy(M) && grounded(M); }
  inline bool floating_ice(int M) { return icy(M) && ocean(M); }
  //! \brief Ice-free cell (grounded or ocean).
  inline bool ice_free(int M) { return !icy(M); }
  inline bool ice_free_ocean(int M) { return ocean(M) && ice_free(M); }
  inline bool ice_free_land(int M) { return grounded(M) && ice_free(M); }
};

class GeometryCalculator
{
public:
  GeometryCalculator(PetscReal seaLevel, const NCConfigVariable &config)
  {
    sea_level = seaLevel;
    alpha = 1 - config.get("ice_density") / config.get("sea_water_density");
    is_dry_simulation = config.get_flag("is_dry_simulation");
  }

  void compute(IceModelVec2S &in_bed, IceModelVec2S &in_thickness,
               IceModelVec2Int &out_mask, IceModelVec2S &out_bed);

  inline void compute(PetscReal bed, PetscReal thickness,
                      int *out_mask, PetscReal *out_surface) {
    const PetscReal  hgrounded = bed + thickness; // FIXME issue #15
    const PetscReal  hfloating = sea_level + alpha*thickness;

    const bool is_floating = hfloating > hgrounded + 1.0,
      ice_free = thickness < 0.01;

    int mask_result;
    PetscReal surface_result;

    if (is_floating && (!is_dry_simulation)) {
      surface_result = hfloating;

      if (ice_free)
        mask_result = MASK_ICE_FREE_OCEAN;
      else
        mask_result = MASK_FLOATING;
    } else {  // Grounded
      surface_result = hgrounded;

      if (ice_free)
        mask_result = MASK_ICE_FREE_BEDROCK;
      else
        mask_result = MASK_GROUNDED;
    }

    if (out_surface != NULL) *out_surface = surface_result;

    if (out_mask != NULL) *out_mask = mask_result;
  }

  inline int mask(PetscReal bed, PetscReal thickness)
  {
    int result;
    compute(bed, thickness, &result, NULL);
    return result;
  }

  inline PetscReal surface(PetscReal bed, PetscReal thickness)
  {
    PetscReal result;
    compute(bed, thickness, NULL, &result);
    return result;
  }

protected:
  PetscReal alpha, sea_level;
  bool is_dry_simulation;
};

class MaskQuery : private Mask
{
public:
  MaskQuery(IceModelVec2Int &m) : mask(m) {}
  
  inline bool ocean(int i, int j) { return Mask::ocean(mask.as_int(i, j)); }

  inline bool grounded(int i, int j) { return !ocean(i, j); }

  inline bool icy(int i, int j) { return Mask::icy(mask.as_int(i, j)); }

  inline bool grounded_ice(int i, int j) { return Mask::grounded_ice(mask.as_int(i, j)); }

  inline bool floating_ice(int i, int j) { return Mask::floating_ice(mask.as_int(i, j)); }

  inline bool ice_free(int i, int j) { return Mask::ice_free(mask.as_int(i, j)); }

  inline bool ice_free_ocean(int i, int j) { return Mask::ice_free_ocean(mask.as_int(i, j)); }

  inline bool ice_free_land(int i, int j) { return Mask::ice_free_land(mask.as_int(i, j)); }

  //! \brief Ice margin (ice-filled with at least one of four neighbors ice-free).
  inline bool ice_margin(int i, int j)
  {
    return icy(i, j) &&
      (ice_free(i + 1, j) || ice_free(i - 1, j) || ice_free(i, j + 1) || ice_free(i, j - 1));
  }

  //! \brief Ice-free margin (at least one of four neighbors has ice).
  inline bool next_to_ice(int i, int j)
  {
    return (icy(i + 1, j) || icy(i - 1, j) || icy(i, j + 1) || icy(i, j - 1));
  }

  //! \brief Ice-free margin (at least one of four neighbors has floating ice).
  inline bool next_to_floating_ice(int i, int j)
  {
    return (floating_ice(i + 1, j) || floating_ice(i - 1, j) || floating_ice(i, j + 1) || floating_ice(i, j - 1));
  }

  //! \brief Ice-free margin (at least one of four neighbors has grounded ice).
  inline bool next_to_grounded_ice(int i, int j)
  {
    return (grounded_ice(i + 1, j) || grounded_ice(i - 1, j) || grounded_ice(i, j + 1) || grounded_ice(i, j - 1));
  }

 //! \brief belongs to margin of an ice shelf but has ice free land neighbor.
  inline bool floating_ice_next_to_icefree_land(int i, int j)
  {
    return floating_ice(i, j) &&
      (ice_free_land(i + 1, j) || ice_free_land(i - 1, j) || ice_free_land(i, j + 1) || ice_free_land(i, j - 1));
  }


  inline PetscErrorCode fill_where_grounded(IceModelVec2S &result, const PetscScalar fillval) {
    PetscErrorCode ierr;

    IceGrid *grid = result.get_grid();

    ierr = mask.begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access(); CHKERRQ(ierr);
    for (PetscInt i = grid->xs; i < grid->xs + grid->xm; ++i) {
      for (PetscInt j = grid->ys; j < grid->ys + grid->ym; ++j) {
        if (grounded(i,j))
          result(i,j) = fillval;
      }
    }
    ierr = result.end_access(); CHKERRQ(ierr);
    ierr = mask.end_access(); CHKERRQ(ierr);

    return 0;
  }

  inline PetscErrorCode fill_where_floating(IceModelVec2S &result, const PetscScalar fillval) {
    PetscErrorCode ierr;

    IceGrid *grid = result.get_grid();

    ierr = mask.begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access(); CHKERRQ(ierr);
    for (PetscInt i = grid->xs; i < grid->xs + grid->xm; ++i) {
      for (PetscInt j = grid->ys; j < grid->ys + grid->ym; ++j) {
        if (ocean(i,j))
          result(i,j) = fillval;
      }
    }
    ierr = result.end_access(); CHKERRQ(ierr);
    ierr = mask.end_access(); CHKERRQ(ierr);

    return 0;
  }
protected:
  IceModelVec2Int &mask;
};

#endif /* _MASK_H_ */
