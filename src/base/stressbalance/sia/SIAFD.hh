// Copyright (C) 2004--2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SIAFD_H_
#define _SIAFD_H_

#include "SSB_Modifier.hh"      // derivesfrom SSB_Modifier
#include "PISMDiagnostic.hh"    // derives from Diag

namespace pism {

class BedSmoother;

/** Implements the shallow ice approximation stress balance.
 *
 * Inputs:
 *
 * - ice geometry (thickness, bed elevation, surface elevation, cell
 *   type mask)
 * - ice enthalpy
 * - ice age (could be used to compute the grain size)
 * - sliding velocity
 *
 * Outputs:
 *
 * - horizontal velocity (3D fields)
 * - diffusive ice flux (for use in the geometry update)
 * - maximum diffusivity (used to determine the maximum allowed time
 *   step length)
 * - volumetric strain heating
 */
class SIAFD : public SSB_Modifier
{
  friend class SIAFD_schoofs_theta;
  friend class SIAFD_topgsmooth;
  friend class SIAFD_thksmooth;
  friend class SIAFD_diffusivity;
  friend class SIAFD_diffusivity_staggered;
  friend class SIAFD_h_x;
  friend class SIAFD_h_y;
public:
  SIAFD(IceGrid &g, EnthalpyConverter &e);

  virtual ~SIAFD();

  virtual void init(Vars &vars);

  virtual void update(IceModelVec2V *vel_input, bool fast);

  //! Add pointers to diagnostic quantities to a dictionary.
  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict);

  virtual void add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &/*result*/) {
  }

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual void define_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                          IO_Type /*nctype*/) {
  }

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual void write_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/) {
  }

protected:

  virtual void compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual void surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual void surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual void surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual void compute_diffusive_flux(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y,
                                                IceModelVec2Stag &result, bool fast);

  virtual void compute_3d_horizontal_velocity(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y,
                                                        IceModelVec2V *vel_input,
                                                        IceModelVec3 &u_out, IceModelVec3 &v_out);

  virtual void compute_I();

  virtual double grainSizeVostok(double age) const;

  virtual void compute_diffusivity(IceModelVec2S &result);
  virtual void compute_diffusivity_staggered(IceModelVec2Stag &result);

  // pointers to input fields:
  IceModelVec2S *bed, *thickness, *surface;
  IceModelVec2Int *mask;
  IceModelVec3 *age, *enthalpy;

  // temporary storage:
  IceModelVec2S work_2d[2];         // for eta, theta and the smoothed thickness
  IceModelVec2Stag work_2d_stag[2]; // for the surface gradient
  IceModelVec3 delta[2];            // store delta on the staggered grid
  IceModelVec3 work_3d[2];      // used to store I and strain_heating
  // on the staggered grid

  BedSmoother *bed_smoother;
  int bed_state_counter;

  // profiling
  int event_sia;

  // unit conversion
  double second_to_kiloyear;
};

//! \brief Computes the multiplier \f$\theta\f$ in Schoof's (2003) theory of the
//! effect of bed roughness on the diffusivity of the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
*/
class SIAFD_schoofs_theta : public Diag<SIAFD>
{
public:
  SIAFD_schoofs_theta(SIAFD *m, IceGrid &g, Vars &my_vars);
  virtual void compute(IceModelVec* &result);
};

//! \brief Computes the smoothed bed elevation from Schoof's (2003) theory of the
//! effect of bed roughness on the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
*/
class SIAFD_topgsmooth : public Diag<SIAFD>
{
public:
  SIAFD_topgsmooth(SIAFD *m, IceGrid &g, Vars &my_vars);
  virtual void compute(IceModelVec* &result);
};

//! \brief Computes the thickness relative to the smoothed bed elevation in
//! Schoof's (2003) theory of the effect of bed roughness on the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
*/
class SIAFD_thksmooth : public Diag<SIAFD>
{
public:
  SIAFD_thksmooth(SIAFD *m, IceGrid &g, Vars &my_vars);
  virtual void compute(IceModelVec* &result);
};

//! \brief Compute diffusivity of the SIA flow.
class SIAFD_diffusivity : public Diag<SIAFD>
{
public:
  SIAFD_diffusivity(SIAFD *m, IceGrid &g, Vars &my_vars);
  virtual void compute(IceModelVec* &result);
};

//! \brief Compute diffusivity of the SIA flow (on the staggered grid).
class SIAFD_diffusivity_staggered : public Diag<SIAFD>
{
public:
  SIAFD_diffusivity_staggered(SIAFD *m, IceGrid &g, Vars &my_vars);
  virtual void compute(IceModelVec* &result);
};

//! \brief Reports the x-component of the ice surface gradient on the staggered
//! grid as computed by SIAFD.
class SIAFD_h_x : public Diag<SIAFD>
{
public:
  SIAFD_h_x(SIAFD *m, IceGrid &g, Vars &my_vars);
  virtual void compute(IceModelVec* &result);
};

//! \brief Reports the y-component of the ice surface gradient on the staggered
//! grid as computed by SIAFD.
class SIAFD_h_y : public Diag<SIAFD>
{
public:
  SIAFD_h_y(SIAFD *m, IceGrid &g, Vars &my_vars);
  virtual void compute(IceModelVec* &result);
};

} // end of namespace pism

#endif /* _SIAFD_H_ */
