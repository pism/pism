// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Constantine Khroulev
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

#include "VariableMetadata.hh"
#include "Timeseries.hh"        // inline code
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "PISMConfigInterface.hh"
#include "iceModelVec.hh"

namespace pism {

//! @brief Class representing diagnostic computations in PISM.
/*!
 * The main goal of this abstraction is to allow accessing metadata
 * corresponding to a diagnostic quantity *before* it is computed.
 *
 * Another goal is to create an interface for computing diagnostics *without*
 * knowing which PISM module is responsible for the computation.
 *
 * Technical note: to compute some diagnostic quantities we need access to
 * protected members of classes. C++ forbids obtaining pointers to non-static
 * methods of a class, but it is possible to define a (friend) function
 *
 * @code
 * IceModelVec::Ptr compute_bar(Foo* model, ...);
 * @endcode
 *
 * which is the same as creating a method `Foo::compute_bar()`, but you *can*
 * get a pointer to it.
 *
 * Diagnostic creates a common interface for all these compute_bar
 * functions.
 */
class Diagnostic {
public:
  Diagnostic(const IceGrid &g);
  virtual ~Diagnostic();

  virtual void update_cumulative();

  //! @brief Compute a diagnostic quantity and return a pointer to a newly-allocated
  //! IceModelVec.
  virtual IceModelVec::Ptr compute() = 0;

  virtual int get_nvars();

  virtual void set_zlevels(std::vector<double> &zlevels);

  virtual SpatialVariableMetadata get_metadata(int N = 0);

  virtual void define(const PIO &nc);

  void set_attrs(const std::string &my_long_name,
                 const std::string &my_standard_name,
                 const std::string &my_units,
                 const std::string &my_glaciological_units,
                 int N = 0);
protected:
  const IceGrid &m_grid;                //!< the grid
  int m_dof;                      //!< number of degrees of freedom; 1 for scalar fields, 2 for vector fields
  IO_Type m_output_datatype;      //!< data type to use in the file
  std::vector<SpatialVariableMetadata> m_vars; //!< metadata corresponding to NetCDF variables
};

//! A template derived from Diagnostic, adding a "Model".
template <class Model>
class Diag : public Diagnostic {
public:
  Diag(Model *m)
    : Diagnostic(m->grid()), model(m) {}
protected:
  Model *model;
};

//! @brief PISM's scalar time-series diagnostics.
class TSDiagnostic {
public:
  TSDiagnostic(const IceGrid &g)
    : m_grid(g), m_ts(NULL) {
  }

  virtual ~TSDiagnostic() {
    delete m_ts;
  }

  virtual void update(double a, double b) = 0;

  virtual void save(double a, double b) {
    if (m_ts) {
      m_ts->interp(a, b);
    }
  }

  virtual void flush() {
    if (m_ts) {
      m_ts->flush();
    }
  }

  virtual void init(const std::string &filename) {
    if (m_ts) {
      m_ts->init(filename);
    }
  }

  virtual std::string get_string(const std::string &name) {
    return m_ts->metadata().get_string(name);
  }

protected:
  const IceGrid &m_grid;                //!< the grid
  DiagnosticTimeseries *m_ts;
};

template <class Model>
class TSDiag : public TSDiagnostic {
public:
  TSDiag(Model *m)
    : TSDiagnostic(m->grid()), model(m) {
    m_time_units = m_grid.ctx()->time()->CF_units_string();
    m_time_dimension_name = m_grid.config()->get_string("time_dimension_name");
  }
protected:
  Model *model;
  std::string m_time_units, m_time_dimension_name;
};

} // end of namespace pism

#endif /* __PISMDiagnostic_hh */
