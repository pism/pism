// Copyright (C) 2011, 2012 PISM Authors
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

#ifndef _PISMMOHRCOULOMBYIELDSTRESS_H_
#define _PISMMOHRCOULOMBYIELDSTRESS_H_

#include "PISMDiagnostic.hh"
#include "PISMYieldStress.hh"
#include "PISMHydrology.hh"
#include "iceModelVec.hh"


//! \brief PISM's default basal yield stress model which applies the Mohr-Coulomb model of deformable, pressurized till.
class PISMMohrCoulombYieldStress : public PISMYieldStress
{
public:
  PISMMohrCoulombYieldStress(IceGrid &g, const NCConfigVariable &conf, PISMHydrology *hydro)
    : PISMYieldStress(g, conf)
  {
    bed_topography = NULL;
    mask = NULL;

    hydrology = hydro;

    if (allocate() != 0) {
      PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMYieldStress constructor.\n");
      PISMEnd();
    }

    ice_density = config.get("ice_density");
    standard_gravity = config.get("standard_gravity");
    till_c_0 = config.get("till_c_0", "kPa", "Pa");
  }

  virtual ~PISMMohrCoulombYieldStress() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);

  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);

  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result);

protected:
  PetscReal standard_gravity, ice_density, till_c_0;
  IceModelVec2S till_phi, tauc, bwp, Po;
  IceModelVec2S *bed_topography;
  IceModelVec2Int *mask;
  PISMVars *variables;
  PISMHydrology *hydrology;

  virtual PetscErrorCode allocate();
  virtual PetscErrorCode topg_to_phi();
  virtual PetscErrorCode tauc_to_phi();
  virtual PetscErrorCode regrid();
};

#endif /* _PISMMOHRCOULOMBYIELDSTRESS_H_ */
