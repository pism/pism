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

#ifndef _PODSLFORCING_H_
#define _PODSLFORCING_H_

#include "PScalarForcing.hh"
#include "POModifier.hh"

namespace pism {

class PO_delta_SL : public PScalarForcing<OceanModel,POModifier>
{
public:
  PO_delta_SL(IceGrid &g, OceanModel* in);
  virtual ~PO_delta_SL();

  virtual void init();
  virtual void sea_level_elevation(double &result);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
protected:
  NCSpatialVariable shelfbmassflux, shelfbtemp;
};

} // end of namespace pism

#endif /* _PODSLFORCING_H_ */
