// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
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

#include <gsl/gsl_math.h>       // GSL_NAN

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
  Diagnostic(IceGrid::ConstPtr g);
  virtual ~Diagnostic();

  typedef std::shared_ptr<Diagnostic> Ptr;

  virtual void update_cumulative();

  //! @brief Compute a diagnostic quantity and return a pointer to a newly-allocated
  //! IceModelVec.
  IceModelVec::Ptr compute();

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
  virtual IceModelVec::Ptr compute_impl() = 0;

  //! the grid
  IceGrid::ConstPtr m_grid;
  //! the unit system
  const units::System::Ptr m_sys;
  //! Configuration flags and parameters
  const Config::ConstPtr m_config;
  //! number of degrees of freedom; 1 for scalar fields, 2 for vector fields
  int m_dof;
  //! data type to use in the file
  IO_Type m_output_datatype;
  //! metadata corresponding to NetCDF variables
  std::vector<SpatialVariableMetadata> m_vars;
  //! fill value (used often enough to justify storing it)
  double m_fill_value;
};

//! A template derived from Diagnostic, adding a "Model".
template <class Model>
class Diag : public Diagnostic {
public:
  Diag(const Model *m)
    : Diagnostic(m->grid()), model(m) {}
protected:
  const Model *model;
};

/*! @brief Helper class for computing time averages of 2D quantities. */
template <class M>
class Diag_average : public Diag<M>
{
public:
  Diag_average(const M *m)
    : Diag<M>(m) {
    m_last_value.create(Diagnostic::m_grid, "last_value_of_cumulative_quantity", WITHOUT_GHOSTS);
    m_last_report_time = GSL_NAN;
  }
protected:
  IceModelVec::Ptr compute_impl() {

    IceModelVec2S::Ptr result(new IceModelVec2S);
    IceModelVec2S &output = *result;
    result->create(Diagnostic::m_grid, "rate_average", WITHOUT_GHOSTS);
    result->metadata(0) = Diagnostic::m_vars[0];

    const IceModelVec2S &value = cumulative_value();
    const double current_time = Diagnostic::m_grid->ctx()->time()->current();

    if (gsl_isnan(m_last_report_time)) {
      const double fill_value = units::convert(Diagnostic::m_sys,
                                               Diagnostic::m_fill_value,
                                               "year-1", "second-1");
      result->set(fill_value);
    } else {
      IceModelVec::AccessList list;
      list.add(*result);
      list.add(m_last_value);
      list.add(value);

      const double dt = current_time - m_last_report_time;

      for (Points p(*Diagnostic::m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        output(i, j) = (value(i, j) - m_last_value(i, j)) / dt;
      }
    }

    // Save the cumulative accumulation and the corresponding time:
    m_last_value.copy_from(value);
    m_last_report_time = current_time;

    return result;
  }

  virtual const IceModelVec2S &cumulative_value() const = 0;

  IceModelVec2S m_last_value;
  double m_last_report_time;
};

//! @brief PISM's scalar time-series diagnostics.
class TSDiagnostic {
public:
  typedef std::shared_ptr<TSDiagnostic> Ptr;

  TSDiagnostic(IceGrid::ConstPtr g)
    : m_grid(g),
      m_config(g->ctx()->config()),
      m_sys(g->ctx()->unit_system()), m_ts(NULL) {
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
  //! the grid
  IceGrid::ConstPtr m_grid;
  //! Configuration flags and parameters
  const Config::ConstPtr m_config;
  //! the unit system
  const units::System::Ptr m_sys;
  DiagnosticTimeseries *m_ts;
};

template <class Model>
class TSDiag : public TSDiagnostic {
public:
  TSDiag(const Model *m)
    : TSDiagnostic(m->grid()), model(m) {
    m_time_units = m_grid->ctx()->time()->CF_units_string();
    m_time_dimension_name = m_grid->ctx()->config()->get_string("time.dimension_name");
  }
protected:
  const Model *model;
  std::string m_time_units, m_time_dimension_name;
};

} // end of namespace pism

#endif /* __PISMDiagnostic_hh */
