// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#ifndef _PISMCONSTANTYIELDSTRESS_H_
#define _PISMCONSTANTYIELDSTRESS_H_

#include "PISMYieldStress.hh"
#include "iceModelVec.hh"
#include "IceGrid.hh"

class PISMConstantYieldStress : public PISMYieldStress
{
public:
  PISMConstantYieldStress(IceGrid &g, const PISMConfig &conf)
    : PISMYieldStress(g, conf)
  {
    if (allocate() != 0) {
      PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMConstantYieldStress constructor.\n");
      PISMEnd();
    }
  }
  virtual ~PISMConstantYieldStress() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result);

  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);

  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result);
protected:
  IceModelVec2S tauc;
  virtual PetscErrorCode allocate();
};

#endif /* _PISMCONSTANTYIELDSTRESS_H_ */
