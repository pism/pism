// Copyright (C) 2004--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _PISMYIELDSTRESS_H_
#define _PISMYIELDSTRESS_H_

#include "PISMComponent.hh"
#include "iceModelVec.hh"

//! Local copy of parameters used by IceModel::getBasalWaterPressure().
struct BWPparams {
  bool usebmr,
    usethkeff;
  PetscReal bmr_scale,
    thkeff_reduce,
    thkeff_H_high,
    thkeff_H_low,
    ice_density,
    standard_gravity;
};

//! \brief The PISM basal yield stress model interface (virtual base class)
class PISMYieldStress : public PISMComponent_TS
{
public:
  PISMYieldStress(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf) {}

  virtual ~PISMYieldStress() {}

  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result) = 0;
};

//! \brief PISM's default basal yield stress model.
class PISMDefaultYieldStress : public PISMYieldStress
{
  friend class PYS_bwp;
public:
  PISMDefaultYieldStress(IceGrid &g, const NCConfigVariable &conf)
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
    p.ice_density   = config.get("ice_density");
    p.standard_gravity = config.get("standard_gravity");
  }

  virtual ~PISMDefaultYieldStress() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, set<string> &result);

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict);

  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc,
                                          nc_type nctype);

  virtual PetscErrorCode write_variables(set<string> vars, string filename);

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result);
protected:
  PetscReal sliding_scale;
  IceModelVec2S till_phi;
  IceModelVec2S *basal_water_thickness, *basal_melt_rate, *ice_thickness,
    *bed_topography;
  IceModelVec2Int *mask;
  BWPparams p;
  PISMVars *variables;

  virtual PetscErrorCode allocate();
  virtual PetscErrorCode topg_to_phi();
  virtual PetscErrorCode regrid();
  virtual PetscScalar basal_water_pressure(PetscScalar thk, PetscScalar bwat,
                                           PetscScalar bmr, PetscScalar frac,
                                           PetscScalar bwat_max);
};

class PISMConstantYieldStress : public PISMYieldStress
{
public:
  PISMConstantYieldStress(IceGrid &g, const NCConfigVariable &conf)
    : PISMYieldStress(g, conf)
  {
    if (allocate() != 0) {
      PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMConstantYieldStress constructor.\n");
      PISMEnd();
    }
  }
  virtual ~PISMConstantYieldStress() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, set<string> &result);

  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc,
                                          nc_type nctype);

  virtual PetscErrorCode write_variables(set<string> vars, string filename);

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result);
protected:
  IceModelVec2S tauc;
  virtual PetscErrorCode allocate();
  virtual PetscErrorCode regrid();
};

//! \brief Computes basal (pore) water pressure using a highly-simplified model.
class PYS_bwp : public PISMDiag<PISMDefaultYieldStress>
{
public:
  PYS_bwp(PISMDefaultYieldStress *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

#endif /* _PISMYIELDSTRESS_H_ */
