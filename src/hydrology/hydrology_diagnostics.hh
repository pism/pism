// Copyright (C) 2012-2017 PISM Authors
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

#include "pism/util/iceModelVec.hh"
#include "pism/util/Diagnostic.hh"
#include "Hydrology.hh"

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
class BasalWaterThickness : public Diag<Routing>
{
public:
  BasalWaterThickness(const Routing *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};


//! \brief Reports the pressure of the transportable water in the subglacial layer.
class BasalWaterPressure : public Diag<Routing>
{
public:
  BasalWaterPressure(const Routing *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};


//! \brief Reports the pressure of the transportable water in the subglacial layer as a fraction of the overburden pressure.
class RelativeBasalWaterPressure : public Diag<Routing>
{
public:
  RelativeBasalWaterPressure(const Routing *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};


//! \brief Reports the effective pressure of the transportable water in the subglacial
//! layer, that is, the overburden pressure minus the pressure.
class EffectiveBasalWaterPressure : public Diag<Routing>
{
public:
  EffectiveBasalWaterPressure(const Routing *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};


//! \brief Report the wall melt rate from dissipation of the potential energy of the
//! transportable water.
class WallMelt : public Diag<Routing>
{
public:
  WallMelt(const Routing *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};


//! \brief Diagnostically reports the staggered-grid components of the velocity of the water in the subglacial layer.
/*! Only available for hydrology::Routing and its derived classes. */
class BasalWaterVelocity : public Diag<Routing>
{
public:
  BasalWaterVelocity(const Routing *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Reports the values of velbase_mag seen by the Hydrology model.
/*! Only available for hydrology::Distributed. */
class Distributed_hydrovelbase_mag : public Diag<Distributed>
{
public:
  Distributed_hydrovelbase_mag(const Distributed *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};


// Diagnostic time-series for mass-conserving ("MC") subglacial hydrology models.
// These eight report the quantities computed in hydrology::Routing::boundary_mass_changes()

//! \brief Reports the rate of loss of liquid water, in kg/s, to locations with mask "ice_free_land()==true".
class MCHydrology_ice_free_land_loss : public TSDiag<TSFluxDiagnostic, Routing>
{
public:
  MCHydrology_ice_free_land_loss(const Routing *m);
  double compute();
};

//! \brief Reports the rate of loss of liquid water, in kg/s, to locations with mask "ocean()==true".
class MCHydrology_ocean_loss : public TSDiag<TSFluxDiagnostic, Routing>
{
public:
  MCHydrology_ocean_loss(const Routing *m);
  double compute();
};

//! \brief Reports the rate of non-conserving gain of liquid water, in kg/s, from water thickness coming out negative during a time step, and being projected up to zero.
class MCHydrology_negative_thickness_gain : public TSDiag<TSFluxDiagnostic, Routing>
{
public:
  MCHydrology_negative_thickness_gain(const Routing *m);
  double compute();
};

//! \brief Reports the rate of loss of liquid water, in kg/s, to locations in the null strip, if that strip has positive width.
class MCHydrology_null_strip_loss : public TSDiag<TSFluxDiagnostic, Routing>
{
public:
  MCHydrology_null_strip_loss(const Routing *m);
  double compute();
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* _HYDROLOGY_DIAGNOSTICS_H_ */

