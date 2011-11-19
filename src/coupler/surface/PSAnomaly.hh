// Copyright (C) 2011 PISM Authors
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

#ifndef _PSANOMALY_H_
#define _PSANOMALY_H_

#include "PGivenClimate.hh"
#include "PSModifier.hh"

//! \brief Reads and uses acab and artm \b anomalies from a file.
class PSAnomaly : public PGivenClimate<PSModifier,PISMSurfaceModel>
{
public:
  PSAnomaly(IceGrid &g, const NCConfigVariable &conf, PISMSurfaceModel* in)
    : PGivenClimate<PSModifier,PISMSurfaceModel>(g, conf, in)
  {
    temp_name = "artm_anomaly";
    mass_flux_name  = "acab_anomaly";
    option_prefix = "-surface_anomaly";
  }
  virtual ~PSAnomaly() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
};

#endif /* _PSANOMALY_H_ */
