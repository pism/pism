// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev
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

#ifndef __PISMDiagnostic_hh
#define __PISMDiagnostic_hh

#include "NCVariable.hh"
#include "Timeseries.hh"        // inline code
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "PISMConfig.hh"
#include "error_handling.hh"

namespace pism {

class IceModelVec;
class Vars;

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
 * Diagnostic creates a common interface for all these compute_bar
 * functions.
 */
class Diagnostic
{
public:
  Diagnostic(IceGrid &g, Vars &my_vars)
    : variables(my_vars), grid(g) {
    output_datatype = PISM_FLOAT;
    dof = 1;
  }
  virtual ~Diagnostic() {}

  //! \brief Update a cumulative quantity needed to compute a rate of change.
  //! So far we there is only one such quantity: the rate of change of the ice
  //! thickness.
  virtual void update_cumulative()
  {
  }

  //! \brief Compute a diagnostic quantity and return a pointer to a newly-allocated
  //! IceModelVec. NB: The caller needs to de-allocate it.
  virtual void compute(IceModelVec* &result) = 0;

  //! Get the number of NetCDF variables corresponding to a diagnostic quantity.
  virtual int get_nvars() {
    return dof;
  }

  //! Reset vertical levels corresponding to the z dimension.
  /** This is called after the automatic grid extension.
   */
  virtual void set_zlevels(std::vector<double> &zlevels)
  {
    for (int j = 0; j < dof; ++j) {
      if (vars[j].get_z().get_name() == "z") {
        vars[j].set_levels(zlevels);
      }
    }
  }

  //! Get a pointer to a metadata object corresponding to variable number N.
  virtual NCSpatialVariable get_metadata(int N = 0)
  {
    if (N >= dof) {
      return NCSpatialVariable(grid.config.get_unit_system(), "missing", grid);
    }

    return vars[N];
  }

  //! Define NetCDF variables corresponding to a diagnostic quantity.
  virtual void define(const PIO &nc)
  {
    for (int j = 0; j < dof; ++j) {
      vars[j].define(nc, output_datatype, true);
    }
  }

  //! \brief A method for setting common variable attributes.
  PetscErrorCode set_attrs(std::string my_long_name,
                           std::string my_standard_name,
                           std::string my_units,
                           std::string my_glaciological_units,
                           int N = 0) {
    PetscErrorCode ierr;

    if (N >= dof) {
      throw RuntimeError::formatted("invalid N (>= dof)");
    }

    vars[N].set_string("pism_intent", "diagnostic");

    vars[N].set_string("long_name", my_long_name);

    vars[N].set_string("standard_name", my_standard_name);

    ierr = vars[N].set_units(my_units); CHKERRQ(ierr);

    if (my_glaciological_units != "") {
      ierr = vars[N].set_glaciological_units(my_glaciological_units); CHKERRQ(ierr);
    }

    return 0;
  }
protected:
  Vars &variables;          //!< dictionary of variables
  IceGrid &grid;                //!< the grid
  int dof;                      //!< number of degrees of freedom; 1 for scalar fields, 2 for vector fields
  IO_Type output_datatype;      //!< data type to use in the file
  std::vector<NCSpatialVariable> vars; //!< metadata corresponding to NetCDF variables
};

//! A template derived from Diagnostic, adding a "Model".
template <class Model>
class Diag : public Diagnostic
{
public:
  Diag(Model *m, IceGrid &g, Vars &my_vars)
    : Diagnostic(g, my_vars), model(m) {}
protected:
  Model *model;
};

//! \brief PISM's scalar time-series diagnostics.
class TSDiagnostic
{
public:
  TSDiagnostic(IceGrid &g, Vars &my_vars)
    : variables(my_vars), grid(g), ts(NULL) {
  }

  virtual ~TSDiagnostic() {
    delete ts;
  }

  virtual void update(double a, double b) = 0;

  virtual void save(double a, double b) {
    if (ts) {
      ts->interp(a, b);
    }
  }

  virtual void flush() {
    if (ts) {
      ts->flush();
    }
  }

  virtual void init(std::string filename) {
    if (ts) {
      ts->init(filename);
    }
  }

  virtual std::string get_string(std::string name) {
    return ts->get_metadata().get_string(name);
  }

protected:
  Vars &variables;          //!< dictionary of variables
  IceGrid &grid;                //!< the grid
  DiagnosticTimeseries *ts;
};

template <class Model>
class TSDiag : public TSDiagnostic
{
public:
  TSDiag(Model *m, IceGrid &g, Vars &my_vars)
    : TSDiagnostic(g, my_vars), model(m) {
    time_units = grid.time->CF_units_string();
    time_dimension_name = grid.config.get_string("time_dimension_name");
  }
protected:
  Model *model;
  std::string time_units, time_dimension_name;
};

} // end of namespace pism

#endif /* __PISMDiagnostic_hh */
