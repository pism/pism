// Copyright (C) 2011 Ed Bueler
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

#ifndef _PISMBEDTHERMALUNIT_H_
#define _PISMBEDTHERMALUNIT_H_

#include "PISMComponent.hh"
#include "iceModelVec.hh"
#include "PISMVars.hh"
#include "materials.hh"
#include "enthalpyConverter.hh"
#include "PISMDiagnostic.hh"

//! Class for a 3d DA-based Vec for PISMBedThermalUnit.
class IceModelVec3BTU : public IceModelVec3D {
public:
  IceModelVec3BTU() { Mbz = -1; Lbz = -1.0; }

  virtual ~IceModelVec3BTU() {}

  virtual PetscErrorCode create(IceGrid &mygrid, const char my_short_name[], bool local,
                                int stencil_width = 1);
                                
  virtual PetscErrorCode get_levels(PetscInt &levels);
  virtual PetscErrorCode get_spacing(PetscReal &dzb);

  PetscErrorCode  stopIfNotLegalLevel(PetscScalar z);

private:
  PetscInt  Mbz;
  PetscReal Lbz;
};


//! Given the ice/bedrock interface temperature for the duration of one time-step, provides upward geothermal flux at that interface.
/*!
The geothermal flux actually applied to the base of an ice sheet is dependent, over time,
on the temperature of the basal ice itself.  The purpose of a bedrock thermal layer
in an ice sheet model is to implement this dependency by using a physical model
for the temperature within that layer, the upper lithosphere.  Because the
upper part of the lithosphere stores or releases energy into the ice,
the typical lithosphere geothermal flux rate is not the same thing as the
geothermal flux applied to the base of the ice.  We regard the lithosphere
geothermal flux rate, which is used in this model as the base of the bedrock
thermal layer, as a time-independent quantity, as it is in all published ice
sheet models.

Let \f$G\f$ be the lithosphere geothermal flux rate, namely
the PISM input variable \c bheatflx; see Related Page \ref std_names .  Let
\f$T_b(t,x,y,z)\f$ be the temperature of the bedrock layer, for elevations
\f$-L_b \le z \le 0\f$; \f$z=0\f$ is the top of the bedrock, the ice/bedrock
interface.  Let \f$k_b\f$ (= \c bedrock_thermal_conductivity in pism_config.cdl)
be the thermal conductivity of the upper lithosphere.  In these terms the actual
upward heat flux into the ice/bedrock interface is the quantity,
  \f[G_0 = -k_b \frac{\partial T_b}{\partial z}.\f]
This is the \e output of the current class, more precisely it is the output of
the method .

 */
class PISMBedThermalUnit : public PISMComponent_TS {

public:
  PISMBedThermalUnit(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &conf);

  virtual ~PISMBedThermalUnit() { }

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);  
  virtual PetscErrorCode write_variables(set<string> vars, string filename);

  virtual PetscErrorCode max_timestep(PetscReal /*t_years*/, PetscReal &dt_years);

  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);

  virtual PetscErrorCode get_upward_geothermal_flux(IceModelVec2S &result);

  IceModelVec3BTU  temp;     //!< storage for bedrock thermal layer temperature;
                             //!    part of state; units K; equally-spaced layers;
                             //!    FIXME: do we want it public?

protected:
  virtual PetscErrorCode allocate();

  IceModelVec2S  ice_base_temp;  //!< temporary storage for boundary value Tb(z=0)

  // parameters of the heat equation:  T_t = D T_xx  where D = k / (rho c)
  PetscScalar    bed_rho, bed_c, bed_k, bed_D;

  EnthalpyConverter &EC; //!< needed to extract base temperature from ice enthalpy

  IceModelVec3   *enthalpy;
  IceModelVec2S  *thk, *ghf;
};

#endif /* _PISMBEDTHERMALUNIT_H_ */

