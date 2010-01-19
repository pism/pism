// Copyright (C) 2009-2010 Andy Aschwanden, Ed Bueler and Constantine Khroulev
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

#ifndef __iceScandModel_hh
#define __iceScandModel_hh

#include <petscda.h>
#include "../../base/grid.hh"
#include "../iceEISModel.hh"

class IceScandModel : public IceEISModel {
public:
    IceScandModel(IceGrid &g, NCConfigVariable &config, NCConfigVariable &overrides);
    virtual PetscErrorCode setFromOptions();
    virtual PetscErrorCode init_couplers();
    
protected:
    virtual PetscErrorCode set_expername_from_options();

    PetscScalar R_cts;
};

#endif /* __iceScandModel_hh */

