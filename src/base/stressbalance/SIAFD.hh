// Copyright (C) 2004--2010 Jed Brown, Ed bueler and Constantine Khroulev
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

#ifndef _SIAFD_H_
#define _SIAFD_H_

#include "SSB_Modifier.hh"
#include "PISMBedSmoother.hh"

class SIAFD : public SSB_Modifier
{
public:
  SIAFD();
  virtual ~SIAFD();

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(IceModelVec2V &vel_input,
                                IceModelVec2S &D2_input,
                                bool fast);

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(PetscInt old_Mz);
protected:
  virtual PetscErrorCode compute_sigma(IceModelVec3 &Sigma);

  virtual PetscErrorCode compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual PetscErrorCode surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual PetscErrorCode surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual PetscErrorCode surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual PetscErrorCode compute_diffusive_flux(IceModelVec2Stag &result, bool fast);

  virtual PetscErrorCode compute_3d_horizontal_velocity(IceModelVec3 &u_out, IceModelVec3 &v_out);

  PetscScalar grainSizeVostok(PetscScalar age) const;

  IceModelVec2S *bed, *thickness, *surface, tmp1, tmp2;
  IceModelVec2Stag tmp3, tmp4;
  IceModelVec3 delta_staggered[2], I_staggered[2], *age, *enthalpy;  

  PISMBedSmoother *sia_bed_smoother;
};


#endif /* _SIAFD_H_ */
