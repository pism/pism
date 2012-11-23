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
#include "iceModelVec.hh"

//! Parameters used by the basal water pressure model.
struct BWPparams {
  bool usebmr,
    usethkeff;
  PetscReal bmr_scale,
    thkeff_reduce,
    thkeff_H_high,
    thkeff_H_low;
};

//! \brief PISM's default basal yield stress model.
class PISMMohrCoulombYieldStress : public PISMYieldStress
{
  friend class PYS_bwp;
public:
  PISMMohrCoulombYieldStress(IceGrid &g, const NCConfigVariable &conf)
    : PISMYieldStress(g, conf)
  {
    sliding_scale = -1.0;
    basal_water_thickness = NULL;
    basal_melt_rate = NULL;
    ice_thickness = NULL;
    bed_topography = NULL;
    mask = NULL;

    if (allocate() != 0) {
      PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMYieldStress constructor.\n");
      PISMEnd();
    }

    p.usebmr        = config.get_flag("bmr_enhance_basal_water_pressure");
    p.usethkeff     = config.get_flag("thk_eff_basal_water_pressure");
    p.bmr_scale     = config.get("bmr_enhance_scale");
    p.thkeff_reduce = config.get("thk_eff_reduced");
    p.thkeff_H_high = config.get("thk_eff_H_high");
    p.thkeff_H_low  = config.get("thk_eff_H_low");

    till_pw_fraction = config.get("till_pw_fraction");
    ice_density = config.get("ice_density");
    standard_gravity = config.get("standard_gravity");
    till_c_0 = config.get("till_c_0", "kPa", "Pa");
    bwat_max = config.get("bwat_max");
  }

  virtual ~PISMMohrCoulombYieldStress() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict);

  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);

  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result);
protected:
  PetscReal standard_gravity, ice_density,
    till_pw_fraction, bwat_max, sliding_scale, till_c_0;
  IceModelVec2S till_phi, tauc;
  IceModelVec2S *basal_water_thickness, *basal_melt_rate, *ice_thickness,
    *bed_topography;
  IceModelVec2Int *mask;
  BWPparams p;
  PISMVars *variables;

  virtual PetscErrorCode allocate();
  virtual PetscErrorCode topg_to_phi();
  virtual PetscErrorCode tauc_to_phi();
  virtual PetscErrorCode regrid();
  virtual PetscReal basal_water_pressure(PetscReal p_overburden, PetscReal bwat,
                                         PetscReal bmr, PetscReal thk);
  virtual PetscReal effective_pressure_on_till(PetscReal p_overburden,
                                               PetscReal b_basal_water);
};

//! \brief Computes basal (pore) water pressure using a highly-simplified model.
class PYS_bwp : public PISMDiag<PISMMohrCoulombYieldStress>
{
public:
  PYS_bwp(PISMMohrCoulombYieldStress *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

#endif /* _PISMMOHRCOULOMBYIELDSTRESS_H_ */
