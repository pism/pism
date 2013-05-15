/* Copyright (C) 2013 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _PISMOCEANKILL_H_
#define _PISMOCEANKILL_H_

#include "PISMComponent.hh"
#include "iceModelVec.hh"

class PISMOceanKill : public PISMComponent {
public:
  PISMOceanKill(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PISMOceanKill();

  virtual PetscErrorCode init(PISMVars &vars);
  PetscErrorCode update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness);

protected:
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, const PIO& nc);

  IceModelVec2Int m_ocean_kill_mask;
};

#endif /* _PISMOCEANKILL_H_ */
