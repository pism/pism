/* Copyright (C) 2022 PISM Authors
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

#ifndef PISM_ICEMODELVEC2S_H
#define PISM_ICEMODELVEC2S_H

#include "pism/util/IceModelVec2.hh"

namespace pism {

/** A class for storing and accessing scalar 2D fields.
    IceModelVec2S is just IceModelVec2 with "dof == 1" */
class IceModelVec2S : public IceModelVec2<double> {
public:
  IceModelVec2S(IceGrid::ConstPtr grid, const std::string &name,
                IceModelVecKind ghostedp, int width = 1);

  typedef std::shared_ptr<IceModelVec2S> Ptr;
  typedef std::shared_ptr<const IceModelVec2S> ConstPtr;
};

std::shared_ptr<IceModelVec2S> duplicate(const IceModelVec2S &source);

// Finite-difference shortcuts. They may be slower than hard-coding FD approximations of x
// and y derivatives. Use with care.
double diff_x(const IceModelVec2S &array, int i, int j);
double diff_y(const IceModelVec2S &array, int i, int j);

// These take grid periodicity into account and use one-sided differences at domain edges.
double diff_x_p(const IceModelVec2S &array, int i, int j);
double diff_y_p(const IceModelVec2S &array, int i, int j);

double sum(const IceModelVec2S &input);
double min(const IceModelVec2S &input);
double max(const IceModelVec2S &input);
double absmax(const IceModelVec2S &input);

void apply_mask(const IceModelVec2S &M, double fill, IceModelVec2S &result);

void compute_magnitude(const IceModelVec2S &v_x,
                       const IceModelVec2S &v_y,
                       IceModelVec2S &result);

} // end of namespace pism

#endif /* PISM_ICEMODELVEC2S_H */
