// Copyright (C) 2004--2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
  SIAFD(IceGrid &g, IceFlowLaw &i, EnthalpyConverter &e, const NCConfigVariable &c)
    : SSB_Modifier(g, i, e, c), WIDE_STENCIL(2) { allocate(); }

  virtual ~SIAFD() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(IceModelVec2V *vel_input,
                                IceModelVec2S *D2_input,
                                bool fast);

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(PetscInt old_Mz);
protected:
  virtual PetscErrorCode allocate();

  virtual PetscErrorCode compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual PetscErrorCode surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual PetscErrorCode surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual PetscErrorCode surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual PetscErrorCode compute_diffusive_flux(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y,
                                                IceModelVec2Stag &result, bool fast);

  virtual PetscErrorCode compute_3d_horizontal_velocity(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y,
                                                        IceModelVec2V *vel_input,
                                                        IceModelVec3 &u_out, IceModelVec3 &v_out);

  virtual PetscErrorCode compute_I();
  virtual PetscErrorCode compute_sigma(IceModelVec2S *D2_input, IceModelVec2Stag &h_x,
                                       IceModelVec2Stag &h_y);

  virtual PetscScalar grainSizeVostok(PetscScalar age) const;

  virtual PetscErrorCode compute_diffusivity(IceModelVec2S &result);

  // pointers to input fields:
  IceModelVec2S *bed, *thickness, *surface;
  IceModelVec2Mask *mask;
  IceModelVec3 *age, *enthalpy;

  // temporary storage:
  IceModelVec2S work_2d[2];         // for eta, theta and the smoothed thickness
  IceModelVec2Stag work_2d_stag[2]; // for the surface gradient
  IceModelVec3 work_3d[4];      // replaces old Sigmastag3 and Istag3; used to
                                // store delta, I and Sigma on the staggered grid

  PISMBedSmoother *bed_smoother;
  const PetscInt WIDE_STENCIL;
  int bed_state_counter;
};


#endif /* _SIAFD_H_ */
