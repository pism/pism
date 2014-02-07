// Copyright (C) 2011, 2014 Constantine Khroulev and David Maxwell
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

#include "Mask.hh"

void GeometryCalculator::compute(IceModelVec2S &bed, IceModelVec2S &thickness,
                    IceModelVec2Int &out_mask, IceModelVec2S &out_surface  )
{
  bed.begin_access();
  thickness.begin_access();
  out_mask.begin_access();
  out_surface.begin_access();
  IceGrid *grid = bed.get_grid();
  
  int GHOSTS = grid->max_stencil_width;
  for (int   i = grid->xs - GHOSTS; i < grid->xs+grid->xm + GHOSTS; ++i) {
    for (int j = grid->ys - GHOSTS; j < grid->ys+grid->ym + GHOSTS; ++j) {
      int mask_value;
      compute(bed(i,j),thickness(i,j),&mask_value,&out_surface(i,j));
      out_mask(i,j) = mask_value;
    } // inner for loop (j)
  } // outer for loop (i)
  
  bed.end_access();
  thickness.end_access();
  out_mask.end_access();
  out_surface.end_access();
}
