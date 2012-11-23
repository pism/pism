// Copyright (C) 2011, 2012 Ed Bueler and Constantine Khroulev
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

class PISMVars;

//! Class for a 3d DA-based Vec for PISMBedThermalUnit.
class IceModelVec3BTU : public IceModelVec3D {
public:
  IceModelVec3BTU() { Lbz = -1.0; }
  virtual ~IceModelVec3BTU() {}

  virtual PetscErrorCode create(IceGrid &mygrid, const char my_short_name[], bool local,
                                int myMbz, PetscReal myLbz, int stencil_width = 1);
                                
  virtual PetscErrorCode get_layer_depth(PetscReal &depth); //!< Return -Lbz value.
  virtual PetscErrorCode get_spacing(PetscReal &dzb);
private:
  PetscReal Lbz;
  bool good_init();
};


//! Given the temperature of the top of the bedrock, for the duration of one time-step, provides upward geothermal flux at that interface at the end of the time-step.
/*!
The geothermal flux actually applied to the base of an ice sheet is dependent, over time,
on the temperature of the basal ice itself.  The purpose of a bedrock thermal layer
in an ice sheet model is to implement this dependency by using a physical model
for the temperature within that layer, the upper lithosphere.  Because the
upper part of the lithosphere stores or releases energy into the ice,
the typical lithosphere geothermal flux rate is not the same thing as the
geothermal flux applied to the base of the ice.  This issue has long been
recognized by ice sheet modelers [%e.g. \ref RitzFabreLetreguilly].

For instance, suppose the ice sheet is in a balanced state in which the geothermal
flux deep in the crust is equal to the heat flux into the ice base.  If the
near-surface ice cools from this state then, because the ice temperature gradient
is now greater in magnitude, between the warm bedrock and the cooler ice, the ice
will for some period receive more than the deep geothermal flux rate. Similarly,
if the ice warms from the balanced state then the temperature difference with
the bedrock has become smaller and the magnitude of the ice basal heat flux will
be less than the deep geothermal rate.

We regard the lithosphere geothermal flux rate, which is applied in this model
to the base of the bedrock thermal layer, as a time-independent quantity.  This
concept is the same as in all published ice sheet models, to our knowledge.

Because the relevant layer of bedrock below an ice sheet is typically shallow,
modeling the bedrock temperature is quite simple.
Let \f$T_b(t,x,y,z)\f$ be the temperature of the bedrock layer, for elevations
\f$-L_b \le z \le 0\f$.  In this routine, \f$z=0\f$ refers to the top of the
bedrock, the ice/bedrock interface.  (Note \f$z=0\f$ is the base of the ice in
IceModel, and thus a different location if ice is floating.)
Let \f$G\f$ be the lithosphere geothermal flux rate, namely the PISM input
variable \c bheatflx; see Related Page \ref std_names .  Let \f$k_b\f$
(= \c bedrock_thermal_conductivity in pism_config.cdl) be the constant thermal
conductivity of the upper lithosphere.  In these terms the actual
upward heat flux into the ice/bedrock interface is the quantity,
  \f[G_0 = -k_b \frac{\partial T_b}{\partial z}.\f]
This is the \e output of the method get_upward_geothermal_flux() in this class.

The evolution equation solved in this class, for which a timestep is done by the
update() method, is the standard 1D heat equation
    \f[\rho_b c_b \frac{\partial T_b}{\partial t} = k_b \frac{\partial^2 T_b}{\partial z^2}\f]
where \f$\rho_b\f$ = \c bedrock_thermal_density and \f$c_b\f$ =
\c bedrock_thermal_specific_heat_capacity in pism_config.cdl.

If \c n_levels>=3 then everything is the general case.  The lithospheric temperature
in \c temp is saved in files as \c litho_temp.  The get_upward_geothermal_flux()
method uses second-order differencing to compute the values of \f$G_0\f$.

If \c n_levels<=1 then this object becomes very simplified: there is no internal
state in IceModelVec3BTU temp.  The update() and allocate() methods are null,
and the get_upward_geothermal_flux() method does nothing other than to copy the
field \f$G\f$ = \c bheatflx into \c result.

If \c n_levels==2 then everything is the general case except that 
get_upward_geothermal_flux() method uses first-order differencing to compute the
values of \f$G_0\f$.
 */
class PISMBedThermalUnit : public PISMComponent_TS {

public:
  PISMBedThermalUnit(IceGrid &g, const NCConfigVariable &conf);

  virtual ~PISMBedThermalUnit() { }

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype);  
  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);

  virtual PetscErrorCode max_timestep(PetscReal /*my_t*/, PetscReal &my_dt, bool &restrict);

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode get_upward_geothermal_flux(IceModelVec2S &result);
protected:
  virtual PetscErrorCode allocate(int Mbz, double Lbz);

  virtual PetscErrorCode bootstrap();
  virtual PetscErrorCode regrid();

  IceModelVec3BTU  temp;     //!< storage for bedrock thermal layer temperature;
                             //!    part of state; units K; equally-spaced layers;
                             //!    This IceModelVec is only created if Mbz > 1.

  // parameters of the heat equation:  T_t = D T_xx  where D = k / (rho c)
  PetscScalar    bed_rho, bed_c, bed_k, bed_D;
  
  PetscInt Mbz;
  PetscReal Lbz;

  IceModelVec2S *bedtoptemp, //!< upper boundary temp, owned by the model to which we are attached
                *ghf; //!< lower boundary flux, owned by the model to which we are attached
};

#endif /* _PISMBEDTHERMALUNIT_H_ */

