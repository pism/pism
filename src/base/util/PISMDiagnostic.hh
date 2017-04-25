// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Constantine Khroulev
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

#include <memory>
#include <deque>

#include "VariableMetadata.hh"
#include "Timeseries.hh"        // inline code
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "PISMConfigInterface.hh"
#include "iceModelVec.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/io/io_helpers.hh"

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

  static Ptr wrap(const IceModelVec2S &input);
  static Ptr wrap(const IceModelVec2V &input);

  void update(double dt);
  void reset();

  //! @brief Compute a diagnostic quantity and return a pointer to a newly-allocated IceModelVec.
  IceModelVec::Ptr compute() const;

  unsigned int n_variables() const;

  SpatialVariableMetadata& metadata(unsigned int N = 0);

  void define(const PIO &file, IO_Type default_type) const;

  void init(const PIO &input, unsigned int time);
  void define_state(const PIO &output) const;
  void write_state(const PIO &output) const;
protected:
  virtual void define_impl(const PIO &file, IO_Type default_type) const;
  virtual void init_impl(const PIO &input, unsigned int time);
  virtual void define_state_impl(const PIO &output) const;
  virtual void write_state_impl(const PIO &output) const;

  void set_attrs(const std::string &my_long_name,
                 const std::string &my_standard_name,
                 const std::string &my_units,
                 const std::string &my_glaciological_units,
                 unsigned int N = 0);

  virtual void update_impl(double dt);
  virtual void reset_impl();

  virtual IceModelVec::Ptr compute_impl() const = 0;

  //! the grid
  IceGrid::ConstPtr m_grid;
  //! the unit system
  const units::System::Ptr m_sys;
  //! Configuration flags and parameters
  const Config::ConstPtr m_config;
  //! number of degrees of freedom; 1 for scalar fields, 2 for vector fields
  unsigned int m_dof;
  //! metadata corresponding to NetCDF variables
  std::vector<SpatialVariableMetadata> m_vars;
  //! fill value (used often enough to justify storing it)
  double m_fill_value;
};

/*!
 * Helper template wrapping quantities with dedicated storage in diagnostic classes.
 *
 * Note: Make sure that that created diagnostics don't outlast fields that they wrap (or you'll hav
 * dangling pointers).
 */
template<class T>
class DiagWithDedicatedStorage : public Diagnostic {
public:
  DiagWithDedicatedStorage(const T &input)
    : Diagnostic(input.get_grid()),
      m_input(input)
  {
    m_dof = input.get_ndof();
    for (unsigned int j = 0; j < m_dof; ++j) {
      m_vars.push_back(input.metadata(j));
    }
  }
protected:
  IceModelVec::Ptr compute_impl() const {
    typename T::Ptr result(new T(m_input.get_grid(), "unnamed", WITHOUT_GHOSTS));
    result->set_name(m_input.get_name());
    for (unsigned int k = 0; k < m_dof; ++k) {
      result->metadata(k) = m_vars[k];
    }

    result->copy_from(m_input);

    return result;
  }
  const T &m_input;
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

/*!
 * Report a time-averaged rate of change of a quantity by accumulating changes over several time
 * steps.
 */
template<class M>
class DiagAverageRate : public Diag<M>
{
public:

  enum InputKind {TOTAL_CHANGE = 0, RATE = 1};

  DiagAverageRate(const M *m, const std::string &name, InputKind kind)
    : Diag<M>(m),
    m_factor(1.0),
    m_input_kind(kind),
    m_accumulator(Diagnostic::m_grid, name + "_accumulator", WITHOUT_GHOSTS),
    m_interval_length(0.0),
    m_time_since_reset(name + "_time_since_reset",
                        Diagnostic::m_config->get_string("time.dimension_name"),
                        Diagnostic::m_sys) {

    m_time_since_reset.set_string("units", "seconds");
    m_time_since_reset.set_string("long_name",
                                  "time since " + m_accumulator.get_name() +
                                  " was reset to 0");

    m_accumulator.metadata().set_string("long_name",
                                        "accumulator for the " + name + " diagnostic");

    m_accumulator.set(0.0);
  }
protected:
  void init_impl(const PIO &input, unsigned int time) {
    m_accumulator.read(input, time);
    {
      std::vector<double> data;
      input.get_1d_var(m_time_since_reset.get_name(),
                       time, 1, // start, count
                       data);
      m_interval_length = data[0];
    }
  }

  void define_state_impl(const PIO &output) const {
    m_accumulator.define(output);
    io::define_timeseries(m_time_since_reset, output, PISM_DOUBLE);
  }

  void write_state_impl(const PIO &output) const {
    m_accumulator.write(output);

    const unsigned int
      time_length = output.inq_dimlen(m_time_since_reset.get_dimension_name()),
      t_start = time_length > 0 ? time_length - 1 : 0;
    io::write_timeseries(output, m_time_since_reset, t_start, m_interval_length, PISM_DOUBLE);
  }

  virtual void update_impl(double dt) {
    // Here the "factor" is used to convert units (from m to kg m-2, for example) and (possibly)
    // integrate over the time integral using the rectangle method.

    double factor = m_factor * (m_input_kind == TOTAL_CHANGE ? 1.0 : dt);

    m_accumulator.add(factor, this->model_input());

    m_interval_length += dt;
  }

  virtual void reset_impl() {
    m_accumulator.set(0.0);
    m_interval_length = 0.0;
  }

  IceModelVec::Ptr compute_impl() const {
    IceModelVec2S::Ptr result(new IceModelVec2S(Diagnostic::m_grid,
                                                "diagnostic", WITHOUT_GHOSTS));
    result->metadata(0) = Diagnostic::m_vars[0];

    if (m_interval_length > 0.0) {
      result->copy_from(m_accumulator);
      result->scale(1.0 / m_interval_length);
    } else {
      std::string
        out = Diagnostic::m_vars[0].get_string("glaciological_units"),
        in  = Diagnostic::m_vars[0].get_string("units");
      const double
        fill = convert(Diagnostic::m_sys, Diagnostic::m_fill_value, out, in);
      result->set(fill);
    }

    return result;
  }
protected:
  // constants initialized in the constructor
  double m_factor;
  InputKind m_input_kind;
  // the state (read from and written to files)
  IceModelVec2S m_accumulator;
  // length of the reporting interval, accumulated along with the cumulative quantity
  double m_interval_length;
  TimeseriesMetadata m_time_since_reset;

  // it should be enough to implement the constructor and this method
  virtual const IceModelVec2S& model_input() {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no default implementation");
  }
};

//! @brief PISM's scalar time-series diagnostics.
class TSDiagnostic {
public:
  typedef std::shared_ptr<TSDiagnostic> Ptr;

  TSDiagnostic(IceGrid::ConstPtr g, const std::string &name);
  virtual ~TSDiagnostic();

  void update(double a, double b);

  void flush();

  void init(const std::string &output_filename,
            std::shared_ptr<std::vector<double>> requested_times);

  const VariableMetadata &metadata() const;

protected:
  virtual double compute(double t0, double t1) = 0;

  void evaluate_regular(double t0, double t1);
  void evaluate_rate(double t0, double t1);

  //! the grid
  IceGrid::ConstPtr m_grid;
  //! Configuration flags and parameters
  const Config::ConstPtr m_config;
  //! the unit system
  const units::System::Ptr m_sys;
  //! time series object used to store computed values and metadata
  std::unique_ptr<Timeseries> m_ts;
  //! requested times
  std::shared_ptr<std::vector<double>> m_times;
  //! times and values, for interpolation
  std::deque<double> m_t, m_v;
  double m_accumulator;
  //! index into m_times
  unsigned int m_current_time;
  //! the name of the file to save to (stored here because it is used by flush(), which is called
  //! from update())
  std::string m_output_filename;
};

template <class Model>
class TSDiag : public TSDiagnostic {
public:
  TSDiag(const Model *m, const std::string &name)
    : TSDiagnostic(m->grid(), name), model(m) {
  }
protected:
  const Model *model;
};

} // end of namespace pism

#endif /* __PISMDiagnostic_hh */
