// Copyright (C) 2004--2012 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "SSB_Modifier.hh"      // derivesfrom SSB_Modifier
#include "PISMDiagnostic.hh"    // derives from PISMDiag

class PISMBedSmoother;

class SIAFD : public SSB_Modifier
{
  friend class SIAFD_schoofs_theta;
  friend class SIAFD_topgsmooth;
  friend class SIAFD_thksmooth;
  friend class SIAFD_diffusivity;
  friend class SIAFD_h_x;
  friend class SIAFD_h_y;
public:
  SIAFD(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &c)
    : SSB_Modifier(g, e, c), WIDE_STENCIL(2) { allocate(); }

  virtual ~SIAFD();

  virtual PetscErrorCode init(PISMVars &vars);

  using PISMComponent_Diag::update;
  virtual PetscErrorCode update(IceModelVec2V *vel_input,
                                IceModelVec2S *D2_input,
                                bool fast);

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(PetscInt old_Mz);

  //! Add pointers to diagnostic quantities to a dictionary.
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict);

  virtual void add_vars_to_output(string /*keyword*/,
                                  map<string,NCSpatialVariable> &/*result*/)
  { }

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode define_variables(set<string> /*vars*/, const PIO &/*nc*/,
                                          PISM_IO_Type /*nctype*/)
  { return 0; }

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode write_variables(set<string> /*vars*/, string /*filename*/)
  { return 0; }

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
  IceModelVec2Int *mask;
  IceModelVec3 *age, *enthalpy;

  // temporary storage:
  IceModelVec2S work_2d[2];         // for eta, theta and the smoothed thickness
  IceModelVec2Stag work_2d_stag[2]; // for the surface gradient
  IceModelVec3 delta[2];            // store delta on the staggered grid
  IceModelVec3 work_3d[2];      // replaces old Sigmastag3 and Istag3; used to
                                // store I and Sigma on the staggered grid

  PISMBedSmoother *bed_smoother;
  const PetscInt WIDE_STENCIL;
  int bed_state_counter;

  // profiling
  int event_sia;

  // unit conversion
  PetscReal second_to_kiloyear;
};

//! \brief Computes the multiplier \f$\theta\f$ in Schoof's (2003) theory of the
//! effect of bed roughness on the diffusivity of the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
 */
class SIAFD_schoofs_theta : public PISMDiag<SIAFD>
{
public:
  SIAFD_schoofs_theta(SIAFD *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the smoothed bed elevation from Schoof's (2003) theory of the
//! effect of bed roughness on the SIA.
/*!
See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
 */
class SIAFD_topgsmooth : public PISMDiag<SIAFD>
{
public:
  SIAFD_topgsmooth(SIAFD *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the thickness relative to the smoothed bed elevation in
//! Schoof's (2003) theory of the effect of bed roughness on the SIA.
/*!
See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
 */
class SIAFD_thksmooth : public PISMDiag<SIAFD>
{
public:
  SIAFD_thksmooth(SIAFD *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Compute diffusivity of the SIA flow.
class SIAFD_diffusivity : public PISMDiag<SIAFD>
{
public:
  SIAFD_diffusivity(SIAFD *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Reports the x-component of the ice surface gradient on the staggered
//! grid as computed by SIAFD.
class SIAFD_h_x : public PISMDiag<SIAFD>
{
public:
  SIAFD_h_x(SIAFD *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Reports the y-component of the ice surface gradient on the staggered
//! grid as computed by SIAFD.
class SIAFD_h_y : public PISMDiag<SIAFD>
{
public:
  SIAFD_h_y(SIAFD *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

#endif /* _SIAFD_H_ */
