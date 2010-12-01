// Copyright (C) 2010 Constantine Khroulev
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

#ifndef __PISMDiagnostic_hh
#define __PISMDiagnostic_hh

#include "iceModelVec.hh"
#include "PISMVars.hh"

class PISMDiagnostic
{
public:
  PISMDiagnostic(IceGrid &g, PISMVars &my_vars)
    : variables(my_vars), grid(g) {
    output_datatype = NC_FLOAT;
    dof = 1;
    vars.resize(dof);
  }
  virtual ~PISMDiagnostic() {}

  //! \brief Compute a diagnostic quantity and return a pointer to a newly-allocated
  //! IceModelVec. NB: The caller needs to de-allocate it.
  virtual PetscErrorCode compute(IceModelVec* &result) = 0;
  
  //! Get the number of NetCDF variables corresponding to a diagnostic quantity.
  virtual int get_nvars() { return dof; }

  //! Get a pointer to a metadata object corresponding to variable number N.
  virtual NCSpatialVariable* get_metadata(int N = 0)
  {
    if (N >= dof) return NULL;

    return &(vars[N]);
  }

  //! Define NetCDF variables corresponding to a diagnostic quantity.
  virtual PetscErrorCode define(const NCTool &nc)
  {
    PetscErrorCode ierr;
    int dummy;

    for (int j = 0; j < dof; ++j) {
      ierr = vars[j].define(nc, dummy, output_datatype, true); CHKERRQ(ierr);
    }

    return 0;
  }

  //! \brief A method for setting common variable attributes.
  PetscErrorCode set_attrs(string my_long_name,
                           string my_standard_name, 
                           string my_units,
                           string my_glaciological_units,
                           int N = 0) {
    PetscErrorCode ierr;

    if (N >= dof) SETERRQ(1, "invalid N (>= dof)");

    vars[N].set_string("pism_intent", "diagnostic");
    vars[N].set_string("long_name", my_long_name);
    vars[N].set_string("standard_name", my_standard_name);
    
    ierr = vars[N].set_units(my_units); CHKERRQ(ierr);
    ierr = vars[N].set_glaciological_units(my_glaciological_units); CHKERRQ(ierr);

    return 0; 
  }
protected:
  PISMVars &variables;          //!< dictionary of variables
  IceGrid &grid;                //!< the grid
  int dof;                      //!< number of degrees of freedom; 1 for scalar fields, 2 for vector fields
  nc_type output_datatype;      //!< data type to use in the file
  vector<NCSpatialVariable> vars; //!< metadata corresponding to NetCDF variables
};

//! A template derived from PISMDiagnostic, adding a "Model".
template <class Model>
class PISMDiag : public PISMDiagnostic
{
public:
  PISMDiag(Model *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiagnostic(g, my_vars), model(m) {}
protected:
  Model *model;
};

#endif /* __PISMDiagnostic_hh */
