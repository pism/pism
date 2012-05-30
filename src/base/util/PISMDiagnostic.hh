// Copyright (C) 2010, 2011, 2012 Constantine Khroulev
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

#include "NCSpatialVariable.hh"
#include "Timeseries.hh"        // inline code
#include "PISMTime.hh"
#include "IceGrid.hh"

class IceModelVec;
class PISMVars;

//! \brief Class representing diagnostic computations in PISM.
/*!
 * The main goal of this abstraction is to allow accessing metadata
 * corresponding to a diagnostic quantity \e before it is computed.
 *
 * Another goal is to create an interface for computing diagnostics \e without
 * knowing which PISM module is responsible for the computation.
 *
 * Technical note: to compute some diagnostic quantities we need access to
 * protected members of classes. C++ forbids obtaining pointers to non-static
 * methods of a class, but it is possible to define a (friend) function
 *
 * \code
 * PetscErrorCode compute_bar(Foo* model, ..., IceModelVec* &result);
 * \endcode
 *
 * which is the same as creating a method Foo::compute_bar(), but you \e can
 * get a pointer to it.
 *
 * PISMDiagnostic creates a common interface for all these compute_bar
 * functions.
 */
class PISMDiagnostic
{
public:
  PISMDiagnostic(IceGrid &g, PISMVars &my_vars)
    : variables(my_vars), grid(g) {
    output_datatype = PISM_FLOAT;
    dof = 1;
    vars.resize(dof);
  }
  virtual ~PISMDiagnostic() {}

  //! \brief Update a cumulative quantity needed to compute a rate of change.
  //! So far we there is only one such quantity: the rate of change of the ice
  //! thickness.
  virtual PetscErrorCode update_cumulative()
  {
    return 0;
  }

  //! \brief Compute a diagnostic quantity and return a pointer to a newly-allocated
  //! IceModelVec. NB: The caller needs to de-allocate it.
  virtual PetscErrorCode compute(IceModelVec* &result) = 0;

  //! Get the number of NetCDF variables corresponding to a diagnostic quantity.
  virtual int get_nvars() { return dof; }

  //! Get a pointer to a metadata object corresponding to variable number N.
  virtual NCSpatialVariable get_metadata(int N = 0)
  {
    if (N >= dof) return NCSpatialVariable();

    return vars[N];
  }

  //! Define NetCDF variables corresponding to a diagnostic quantity.
  virtual PetscErrorCode define(const PIO &nc)
  {
    PetscErrorCode ierr;

    for (int j = 0; j < dof; ++j) {
      ierr = vars[j].define(nc, output_datatype, true); CHKERRQ(ierr);
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

    if (N >= dof) SETERRQ(grid.com, 1, "invalid N (>= dof)");

    vars[N].set_string("pism_intent", "diagnostic");
    if (my_long_name != "")
      vars[N].set_string("long_name", my_long_name);

    if (my_standard_name != "")
      vars[N].set_string("standard_name", my_standard_name);

    if (my_units != "") {
      ierr = vars[N].set_units(my_units); CHKERRQ(ierr);
    }

    if (my_glaciological_units != "") {
      ierr = vars[N].set_glaciological_units(my_glaciological_units); CHKERRQ(ierr);
    }

    return 0;
  }
protected:
  PISMVars &variables;          //!< dictionary of variables
  IceGrid &grid;                //!< the grid
  int dof;                      //!< number of degrees of freedom; 1 for scalar fields, 2 for vector fields
  PISM_IO_Type output_datatype;      //!< data type to use in the file
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

//! \brief PISM's scalar time-series diagnostics.
class PISMTSDiagnostic
{
public:
  PISMTSDiagnostic(IceGrid &g, PISMVars &my_vars)
    : variables(my_vars), grid(g), ts(NULL) {
  }

  virtual ~PISMTSDiagnostic() {
    delete ts;
  }

  virtual PetscErrorCode update(PetscReal a, PetscReal b) = 0;

  virtual PetscErrorCode save(PetscReal a, PetscReal b) {
    if (ts)
      return ts->interp(a, b);

    return 0;
  }

  virtual PetscErrorCode flush() {
    if (ts)
      return ts->flush();

    return 0;
  }

  virtual PetscErrorCode init(string filename) {
    if (ts)
      return ts->init(filename);
    return 0;
  }

  virtual string get_string(string name) {
    return ts->get_string(name);
  }

protected:
  PISMVars &variables;          //!< dictionary of variables
  IceGrid &grid;                //!< the grid
  DiagnosticTimeseries *ts;
};

template <class Model>
class PISMTSDiag : public PISMTSDiagnostic
{
public:
  PISMTSDiag(Model *m, IceGrid &g, PISMVars &my_vars)
    : PISMTSDiagnostic(g, my_vars), model(m) {
    time_units = grid.time->CF_units();
    time_dimension_name = grid.config.get_string("time_dimension_name");
  }
protected:
  Model *model;
  string time_units, time_dimension_name;
};

#endif /* __PISMDiagnostic_hh */
