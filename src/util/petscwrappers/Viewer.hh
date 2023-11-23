/* Copyright (C) 2014, 2015 PISM Authors
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

#ifndef _VIEWER_H_
#define _VIEWER_H_

#include <petscviewer.h>
#include <string>

#include "pism/util/Wrapper.hh"

namespace pism {
namespace petsc {
class Viewer : public Wrapper<PetscViewer> {
public:
  Viewer(MPI_Comm com, const std::string &name,
         unsigned int target_size, double Lx, double Ly);
  Viewer(PetscViewer v);
  Viewer(MPI_Comm c);
  Viewer();
  ~Viewer();
private:
  void compute_size(unsigned int target_size, double Lx, double Ly,
                    unsigned int &X, unsigned int &Y);
};
} // end of namespace petsc
} // end of namespace pism

#endif /* _VIEWER_H_ */
