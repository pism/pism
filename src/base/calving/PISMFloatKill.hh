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

#ifndef _PISMFLOATKILL_H_
#define _PISMFLOATKILL_H_

#include "PISMComponent.hh"

namespace pism {

class IceModelVec2Int;
class IceModelVec2S;

/*! \brief Calving mechanism removing floating ice. */
class FloatKill : public Component
{
public:
  FloatKill(IceGrid &g, const Config &conf);
  virtual ~FloatKill();

  virtual void init(Vars &vars);
  void update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO& nc);
};


} // end of namespace pism

#endif /* _PISMFLOATKILL_H_ */
