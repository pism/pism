// Copyright (C) 2012-2016 PISM Authors
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

#ifndef _HYDROLOGY_DIAGNOSTICS_H_
#define _HYDROLOGY_DIAGNOSTICS_H_

#include "base/util/iceModelVec.hh"
#include "base/util/PISMDiagnostic.hh"
#include "PISMHydrology.hh"

namespace pism {
namespace hydrology {
/*! \file
  Interfaces for the following diagnostics which are handled by Hydrology
  instances; some of these may be replaced by state variables; listed by short
  name:
  * bwat [replace by state var in hydrology::Routing and hydrology::Distributed]
  * bwp [replace by state var in hydrology::Distributed]
  * bwprel
  * effbwp
  * hydroinput
  * wallmelt
  Interfaces for the following diagnostics which are handled by
  hydrology::Routing instances:
  * bwatvel
  */


//! \brief Reports the thickness of the transportable water in the subglacial layer.
class Hydrology_bwat : public Diag<Hydrology>
{
public:
  Hydrology_bwat(const Hydrology *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Reports the pressure of the transportable water in the subglacial layer.
class Hydrology_bwp : public Diag<Hydrology>
{
public:
  Hydrology_bwp(const Hydrology *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Reports the pressure of the transportable water in the subglacial layer as a fraction of the overburden pressure.
class Hydrology_bwprel : public Diag<Hydrology>
{
public:
  Hydrology_bwprel(const Hydrology *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Reports the effective pressure of the transportable water in the subglacial layer, that is, the overburden pressure minus the pressure.
class Hydrology_effbwp : public Diag<Hydrology>
{
public:
  Hydrology_effbwp(const Hydrology *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Reports the values of bmelt seen by the Hydrology model.
class Hydrology_hydrobmelt : public Diag<Hydrology>
{
public:
  Hydrology_hydrobmelt(const Hydrology *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Reports the total input rate of water into the subglacial layer.
class Hydrology_hydroinput : public Diag<Hydrology>
{
public:
  Hydrology_hydroinput(const Hydrology *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Report the wall melt rate from dissipation of the potential energy of the transportable water.
class Hydrology_wallmelt : public Diag<Hydrology>
{
public:
  Hydrology_wallmelt(const Hydrology *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Diagnostically reports the staggered-grid components of the velocity of the water in the subglacial layer.
/*! Only available for hydrology::Routing and its derived classes. */
class Routing_bwatvel : public Diag<Routing>
{
public:
  Routing_bwatvel(const Routing *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Reports the values of velbase_mag seen by the Hydrology model.
/*! Only available for hydrology::Distributed. */
class Distributed_hydrovelbase_mag : public Diag<Distributed>
{
public:
  Distributed_hydrovelbase_mag(const Distributed *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


// Diagnostic time-series for mass-conserving ("MC") subglacial hydrology models.
// These eight report the quantities computed in hydrology::Routing::boundary_mass_changes()

//! \brief Reports the cumulative loss of liquid water, in kg, to locations with mask "ice_free_land()==true".
class MCHydrology_ice_free_land_loss_cumulative : public TSDiag<Routing>
{
public:
  MCHydrology_ice_free_land_loss_cumulative(const Routing *m);
  virtual void update(double a, double b);
};

//! \brief Reports the rate of loss of liquid water, in kg/s, to locations with mask "ice_free_land()==true".
class MCHydrology_ice_free_land_loss : public TSDiag<Routing>
{
public:
  MCHydrology_ice_free_land_loss(const Routing *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative loss of liquid water, in kg, to locations with mask "ocean()==true".
class MCHydrology_ocean_loss_cumulative : public TSDiag<Routing>
{
public:
  MCHydrology_ocean_loss_cumulative(const Routing *m);
  virtual void update(double a, double b);
};

//! \brief Reports the rate of loss of liquid water, in kg/s, to locations with mask "ocean()==true".
class MCHydrology_ocean_loss : public TSDiag<Routing>
{
public:
  MCHydrology_ocean_loss(const Routing *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative non-conserving gain of liquid water, in kg, from water thickness coming out negative during a time step, and being projected up to zero.
class MCHydrology_negative_thickness_gain_cumulative : public TSDiag<Routing>
{
public:
  MCHydrology_negative_thickness_gain_cumulative(const Routing *m);
  virtual void update(double a, double b);
};

//! \brief Reports the rate of non-conserving gain of liquid water, in kg/s, from water thickness coming out negative during a time step, and being projected up to zero.
class MCHydrology_negative_thickness_gain : public TSDiag<Routing>
{
public:
  MCHydrology_negative_thickness_gain(const Routing *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative loss of liquid water, in kg, to locations in the null strip, if that strip has positive width.
class MCHydrology_null_strip_loss_cumulative : public TSDiag<Routing>
{
public:
  MCHydrology_null_strip_loss_cumulative(const Routing *m);
  virtual void update(double a, double b);
};

//! \brief Reports the rate of loss of liquid water, in kg/s, to locations in the null strip, if that strip has positive width.
class MCHydrology_null_strip_loss : public TSDiag<Routing>
{
public:
  MCHydrology_null_strip_loss(const Routing *m);
  virtual void update(double a, double b);
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* _HYDROLOGY_DIAGNOSTICS_H_ */

