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

#ifndef _POFACTORY_H_
#define _POFACTORY_H_

#include "PCFactory.hh"
#include "POModifier.hh"

class POFactory : public PCFactory<PISMOceanModel,POModifier> {
public:
  POFactory(IceGrid& g, const NCConfigVariable& conf)
    : PCFactory<PISMOceanModel,POModifier>(g, conf)
  {
    add_standard_types();
    option = "ocean";
  }
  virtual ~POFactory() {}
  virtual void add_standard_types();
};

#endif /* _POFACTORY_H_ */
