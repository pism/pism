// Copyright (C) 2011, 2014 PISM Authors
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

#ifndef _PAFACTORY_H_
#define _PAFACTORY_H_

#include "PISMAtmosphere.hh"
#include "PAModifier.hh"
#include "PCFactory.hh"

namespace pism {

class PAFactory : public PCFactory<AtmosphereModel,PAModifier> {
public:
  PAFactory(IceGrid& g, const Config& conf)
    : PCFactory<AtmosphereModel,PAModifier>(g, conf)
  {
    add_standard_types();
    option = "atmosphere";
  }
  virtual ~PAFactory() {}
  virtual void add_standard_types();
};

} // end of namespace pism

#endif /* _PAFACTORY_H_ */
