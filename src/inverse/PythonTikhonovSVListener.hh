// Copyright (C) 2012  David Maxwell
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

#ifndef PYTHONTIKHONOVSVLISTENER_HH_ZILN5E62
#define PYTHONTIKHONOVSVLISTENER_HH_ZILN5E62

#include "iceModelVec.hh"
#include <memory>

class PythonTikhonovSVListener {
public:
  typedef std::tr1::shared_ptr<PythonTikhonovSVListener> Ptr;
  PythonTikhonovSVListener() {}
  virtual ~PythonTikhonovSVListener() {}
  virtual void iteration(PetscInt iter, PetscReal eta,
   PetscReal objectiveValue, PetscReal designValue,
   IceModelVec2S &d, IceModelVec2S &diff_d, IceModelVec2S &grad_d,
   IceModelVec2V &u,  IceModelVec2V &diff_u,  IceModelVec2S &grad_u,
   IceModelVec2S &gradient) { 
     (void) iter; (void) eta; (void) objectiveValue; (void) designValue;
     (void) d; (void) diff_d; (void) grad_d; 
     (void) u; (void) diff_u; (void) grad_u;
     (void) gradient;
  };
};

#endif /* end of include guard: PYTHONTIKHONOVSVLISTENER_HH_ZILN5E62 */
