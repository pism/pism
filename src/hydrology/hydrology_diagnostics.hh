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
#include "Routing.hh"

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

} // end of namespace hydrology
} // end of namespace pism

#endif /* _HYDROLOGY_DIAGNOSTICS_H_ */

