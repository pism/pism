// Copyright (C) 2012-2014 PISM Authors
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

#include "iceModelVec.hh"
#include "PISMDiagnostic.hh"
#include "PISMHydrology.hh"

namespace pism {
/*! \file
  Interfaces for the following diagnostics which are handled by PISMHydrology
  instances; some of these may be replaced by state variables; listed by short
  name:
  * bwat [replace by state var in PISMRoutingHydrology and PISMDistributedHydrology]
  * bwp [replace by state var in PISMDistributedHydrology]
  * bwprel
  * effbwp
  * hydroinput
  * wallmelt
  Interfaces for the following diagnostics which are handled by
  PISMRoutingHydrology instances:
  * bwatvel
  */


//! \brief Reports the thickness of the transportable water in the subglacial layer.
class PISMHydrology_bwat : public PISMDiag<PISMHydrology>
{
public:
  PISMHydrology_bwat(PISMHydrology *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};


//! \brief Reports the pressure of the transportable water in the subglacial layer.
class PISMHydrology_bwp : public PISMDiag<PISMHydrology>
{
public:
  PISMHydrology_bwp(PISMHydrology *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};


//! \brief Reports the pressure of the transportable water in the subglacial layer as a fraction of the overburden pressure.
class PISMHydrology_bwprel : public PISMDiag<PISMHydrology>
{
public:
  PISMHydrology_bwprel(PISMHydrology *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};


//! \brief Reports the effective pressure of the transportable water in the subglacial layer, that is, the overburden pressure minus the pressure.
class PISMHydrology_effbwp : public PISMDiag<PISMHydrology>
{
public:
  PISMHydrology_effbwp(PISMHydrology *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};


//! \brief Reports the total input rate of water into the subglacial layer.
class PISMHydrology_hydroinput : public PISMDiag<PISMHydrology>
{
public:
  PISMHydrology_hydroinput(PISMHydrology *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};


//! \brief Report the wall melt rate from dissipation of the potential energy of the transportable water.
class PISMHydrology_wallmelt : public PISMDiag<PISMHydrology>
{
public:
  PISMHydrology_wallmelt(PISMHydrology *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};


//! \brief Diagnostically reports the staggered-grid components of the velocity of the water in the subglacial layer.
/*! Only available for PISMRoutingHydrology and its derived classes. */
class PISMRoutingHydrology_bwatvel : public PISMDiag<PISMRoutingHydrology>
{
public:
  PISMRoutingHydrology_bwatvel(PISMRoutingHydrology *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

} // end of namespace pism

#endif /* _HYDROLOGY_DIAGNOSTICS_H_ */

