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

#ifndef _SIA_SLIDING_H_
#define _SIA_SLIDING_H_

#include "ShallowStressBalance.hh"

namespace pism {

/*!
 * This class implements an SIA sliding law.
 *
 * It is used by pismv test E \b only, hence the code duplication (the surface
 * gradient code is from SIAFD).
 */
class SIA_Sliding : public ShallowStressBalance
{
public:
  SIA_Sliding(IceGrid &g, EnthalpyConverter &e, const Config &conf);

  virtual ~SIA_Sliding();

  virtual PetscErrorCode init(Vars &vars);

  virtual PetscErrorCode update(bool fast, IceModelVec2S &melange_back_pressure);

  virtual void add_vars_to_output(std::string /*keyword*/, std::set<std::string> &/*result*/)
  { }

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode define_variables(std::set<std::string> /*vars*/, const PIO &/*nc*/,
                                          IO_Type /*nctype*/)
  { return 0; }

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode write_variables(std::set<std::string> /*vars*/, const PIO &/*nc*/)
  { return 0; }

protected:
  virtual PetscErrorCode allocate();

  virtual PetscErrorCode compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual PetscErrorCode surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual PetscErrorCode surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual PetscErrorCode surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual double basalVelocitySIA(double /*x*/, double /*y*/,
                                  double H, double T,
                                  double /*alpha*/, double mu,
                                  double min_T) const;
  IceModelVec2Int *mask;
  IceModelVec2S *thickness, *surface, *bed, work_2d;
  IceModelVec3 *enthalpy;
  IceModelVec2Stag work_2d_stag[2]; // for the surface gradient
  double standard_gravity;

  bool verification_mode;
  std::string eisII_experiment;
};

} // end of namespace pism

#endif /* _SIA_SLIDING_H_ */
