/* Copyright (C) 2013, 2014 PISM Authors
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

#ifndef _PO_DELTA_MBP_H_
#define _PO_DELTA_MBP_H_

#include "PScalarForcing.hh"
#include "POModifier.hh"

namespace pism {

/**
 * Scalar melange back-pressure fraction forcing.
 * 
 */
class PO_delta_MBP : public PScalarForcing<OceanModel,POModifier>
{
public:
  PO_delta_MBP(IceGrid &g, const Config &conf, OceanModel* in);
  virtual ~PO_delta_MBP();

  virtual PetscErrorCode init(Vars &vars);
  virtual PetscErrorCode melange_back_pressure_fraction(IceModelVec2S &result);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual PetscErrorCode define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual PetscErrorCode write_variables(const std::set<std::string> &vars, const PIO &nc);
protected:
  NCSpatialVariable shelfbmassflux, shelfbtemp;
private:
  PetscErrorCode allocate_PO_delta_MBP();
};

} // end of namespace pism

#endif /* _PO_DELTA_MBP_H_ */
