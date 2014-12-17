// Copyright (C) 2011, 2012, 2014 PISM Authors
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

#ifndef _PSTEMPERATUREINDEX_OLD_H_
#define _PSTEMPERATUREINDEX_OLD_H_

#include "PISMSurface.hh"
#include "localMassBalance_old.hh"
#include "NCVariable.hh"

namespace pism {

//! \brief A class implementing a temperature-index (positive degree-day) scheme
//! to compute melt and runoff, and thus surface mass balance, from
//! precipitation and air temperature.
/*! 
  Temperature-index schemes are far from perfect as a way of modeling surface mass
  balance on ice sheets which experience surface melt, but they are known to have
  reasonable data requirements and to do a good job when tuned appropriately
  [\ref Hock05].

  This base class already accesses a fair amount of functionality.  It holds a
  pointer to an instance of the LocalMassBalance class.  This class has method
  LocalMassBalance::getMassFluxFromTemperatureTimeSeries() which uses the
  precipitation during the ice sheet model time step, plus a variable temperature
  over that time step, to compute melt, refreeze, and surface balance.

  This base class reads options <tt>-pdd_factor_snow</tt>, <tt>-pdd_factor_ice</tt>,
  and <tt>-pdd_refreeze</tt> and sets these factors accordingly, in the case where
  the factors are independent of location.  If option <tt>-pdd_fausto</tt> is used
  then an object is called which updates these values based on the location.
*/
class PSTemperatureIndex_Old : public SurfaceModel {
public:
  PSTemperatureIndex_Old(IceGrid &g);
  virtual ~PSTemperatureIndex_Old();
  virtual void update(PetscReal my_t, PetscReal my_dt);
  virtual void init();
  virtual void max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict);
  virtual void ice_surface_mass_flux(IceModelVec2S &result);
  virtual void ice_surface_temperature(IceModelVec2S &result);
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);  
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
protected:
  virtual void update_internal(PetscReal my_t, PetscReal my_dt);
  LocalMassBalance_Old *mbscheme;	      //!< mass balance scheme to use

  FaustoGrevePDDObject_Old *faustogreve;  //!< if not NULL then user wanted fausto PDD stuff

  DegreeDayFactors_Old base_ddf;          //!< holds degree-day factors in location-independent case
  PetscScalar  base_pddStdDev,        //!< K; daily amount of randomness
    base_pddThresholdTemp; //!< K; temps are positive above this
  IceModelVec2S
  climatic_mass_balance,      //!< cached surface mass balance rate
    accumulation_rate, //!< diagnostic output accumulation rate (snow - rain)
    melt_rate,    //!< diagnostic output melt rate (rate at which snow
  //! and ice is melted, but some snow melt refreezes)
    runoff_rate;          //!< diagnostic output meltwater runoff rate

  IceModelVec2S *lat, *lon, *usurf;  //!< PSTemperatureIndex_Old must hold these
  //!pointers in order to use object which
  //!needs 3D location to determine degree
  //!day factors.
  bool pdd_annualize;
  PetscReal next_pdd_update;

private:
  std::string temperature_name, mass_balance_name;
  NCSpatialVariable ice_surface_temp;
};

} // end of namespace pism

#endif /* _PSTEMPERATUREINDEX_OLD_H_ */
