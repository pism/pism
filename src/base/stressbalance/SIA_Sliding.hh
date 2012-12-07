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

#ifndef _SIA_SLIDING_H_
#define _SIA_SLIDING_H_

#include "ShallowStressBalance.hh"

/*!
 * This class implements an SIA sliding law.
 *
 * It is used by pismv test E \b only, hence the code duplication (the surface
 * gradient code is from SIAFD).
 */
class SIA_Sliding : public ShallowStressBalance
{
public:
  SIA_Sliding(IceGrid &g, IceBasalResistancePlasticLaw &b,
              EnthalpyConverter &e, const NCConfigVariable &conf)
    : ShallowStressBalance(g, b, e, conf)
  {
    verification_mode = false;
    eisII_experiment = "";
    allocate();
  }

  virtual ~SIA_Sliding()
  {
    if (flow_law != NULL) {
      delete flow_law;
      flow_law = NULL;
    }
  }

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(bool fast);

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
  virtual PetscErrorCode write_variables(set<string> /*vars*/, const PIO &/*nc*/)
  { return 0; }

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict) {
    dict["taud"] = new SSB_taud(this, grid, *variables);
    dict["taud_mag"] = new SSB_taud_mag(this, grid, *variables);
  }

protected:
  virtual PetscErrorCode allocate();

  virtual PetscErrorCode compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual PetscErrorCode surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual PetscErrorCode surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual PetscErrorCode surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual PetscScalar basalVelocitySIA(PetscScalar /*x*/, PetscScalar /*y*/,
                                       PetscScalar H, PetscScalar T,
                                       PetscScalar /*alpha*/, PetscScalar mu,
                                       PetscScalar min_T) const;
  IceModelVec2Int *mask;
  IceModelVec2S *thickness, *surface, *bed, work_2d;
  IceModelVec3 *enthalpy;
  IceModelVec2Stag work_2d_stag[2]; // for the surface gradient
  double standard_gravity;

  bool verification_mode;
  string eisII_experiment;
};

#endif /* _SIA_SLIDING_H_ */
