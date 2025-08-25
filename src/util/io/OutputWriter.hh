/* Copyright (C) 2025 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef PISM_OUTPUTWRITER_H
#define PISM_OUTPUTWRITER_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

namespace pism {

class Config;
class SpatialVariableMetadata;
class VariableMetadata;

namespace grid {
class DistributedGridInfo;
}

namespace io {
enum Type : int;
}

namespace units {
class System;
}

/*!
 * File output API
 *
 * PISM writes the following kinds of output:
 *
 * 1) output at the end of a run (used to re-start the model),
 *
 * 2) 2D and 3D diagnostic output saved at specified model times during a run (affects
 *    model time stepping),
 *
 * 3) Same as 2), but each time record is saved to a separate file,
 *
 * 4) Same as 2), but appending to a file created by an earlier run,
 *
 * 5) Snapshots of the model state at times *close to* specified model times during a run
 *    (does not affect model time stepping; can be used to re-start a failed run)
 *
 * 6) Same as 5, but each snapshot is saved to separate file,
 *
 * 7) Scalar time-dependent diagnostics (value are stored *redundantly* on all MPI ranks),
 *
 * 8) Same as 7), but appending to a file created by an earlier run,
 *
 * 9) Snapshots of the model state saved after a specified *wall clock time* interval
 *    passed (used to re-start a failed run).
 *
 * All files contain *one* unlimited dimension (time).
 *
 * File contents are determined at run time.
 *
 * A file may contain a mix of 1D, 2D, and 3D time-dependent and time-independent
 * variables.
 *
 * A file may contain more than one x,y grid and more than one set of vertical (z) levels.
 *
 * Appending to a file requires being able to get the current length of the time dimension
 * in a file and the last value of the corresponding coordinate variable.
 *
 * In this API, the first call using `file_name` opens the file `file_name`. If the file already exists
 * it is moved to `file_name` + "~" (a "backup" file).
 *
 * If the first call using `file_name` is `append(file_name)`, the file is opened for
 * appending.
 *
 * An opened file remains open until `close()` is called or until an instance of
 * `OutputWriter` is de-allocated (i.e. until the end of a model run).
 *
 * PISM defines all variables before writing *any* of the associated data. Attributes are
 * set *once* and not modified afterwards. (This should make it possible to aggregate all
 * metadata and write all of it at once.)
 *
 * PISM buffers scalar time-dependent diagnostics to reduce the number of I/O operations.
 * 2D and 3D arrays are written one time record at a time (increase the length of time
 * dimension by one, write a bunch of variables, increase the length of time dimension by
 * one, write more, etc).
 */
class OutputWriter {
public:
  OutputWriter(MPI_Comm comm, const Config &config);
  virtual ~OutputWriter();

  /*!
   * Add attributes to all 2D and 3D variables using "x" and "y" dimensions except for
   * "lat", "lon", "lat_bnds" and "lon_bnds".
   *
   * These attributes (usually 'coordinates = "lat lon"' and "grid_mapping") are used to
   * provide grid information.
   *
   * A variable should have the "coordinates" attribute only if the file contains
   * variables "lat" and "lon". Similarly, a variable should have the "grid_mapping"
   * attribute only if the file contains a grid mapping variable. The code writing a
   * variable has no way to determine the contents of an output file, so this method
   * allows IceModel to use the information *it* has to adjust metadata of 2D and 3D
   * variables written to output files.
   */
  void add_extra_attributes(const std::string &file_name,
                            const std::map<std::string, std::string> &attributes);

  /*!
   * Define a dimension.
   * 
   * No-op if the dimension already exists.
   *
   * `length` of zero corresponds to an unlimited dimension.
   */
  void define_dimension(const std::string &file_name, const std::string &dimension_name,
                        unsigned int length);

  /*!
   * Define a variable given a list of dimension names and set its attributes.
   *
   * Use this method to define coordinate variables (`x`, `y`, `time`, etc). For scalar
   * time-dependent model outputs, use `define_timeseries_variable()`.
   *
   * No-op if the variable already exists.
   *
   * The name of the variable is obtained using `metadata.get_name()` and its type using `metadata.get_output_type()`.
   */
  void define_variable(const std::string &file_name, const VariableMetadata &metadata,
                       const std::vector<std::string> &dims);

  /*!
   * Add a spatial variable to the list of variables that can be written to output files.
   *
   * This has to be done before a spatial variable is *defined* in an output file.
   *
   * @param[in] metadata variable metadata (name, attributes, etc)
   * @param[in] grid domain decomposition information
   */
  void add_spatial_variable(const SpatialVariableMetadata &metadata,
                            const grid::DistributedGridInfo &grid);

  /*!
   * Define a 2D or 3D (possibly time-dependent) variable and set its attributes.
   *
   * No-op if the variable already exists.
   *
   * Stores domain decomposition `grid` and makes it accessible using grid_info().
   *
   * @param[in] file_name name of the output file
   * @param[in] variable_name variable name
   */
  void define_spatial_variable(const std::string &file_name,
                               const std::string &variable_name);

  /*!
   * Define a scalar time-dependent variable and set its attributes.
   *
   * This method should be used to define variables containing scalar time-dependent model
   * outputs. Use `define_variable()` to define coordinate variables (`x`, `y`, `time`,
   * etc).
   *
   * No-op if the variable already exists.
   *
   * @param[in] file_name name of the output file
   * @param[in] metadata variable metadata (name, attributes, etc)
   */
  void define_timeseries_variable(const std::string &file_name,
                                  const VariableMetadata &metadata);

  /*!
   * Set global attributes for a given output file.
   *
   * Numbers are written as NC_DOUBLE.
   */
  void set_global_attributes(const std::string &file_name,
                             const std::map<std::string, std::string> &strings,
                             const std::map<std::string, std::vector<double> > &numbers);

  /*!
   * Append to the global attribute "history" in the output file.
   *
   * The "history" attribute is treated as a newline-delimited list. This call adds a "\n"
   * followed by the string in `text` to the attribute "history" in the file `file_name`.
   */
  void append_history(const std::string &file_name, const std::string &text);

  /*!
   * Increase the length of the time dimension by one, appending the value `time_seconds`.
   *
   * This should increase the value returned by `time_dimension_length()` by one.
   */
  void append_time(const std::string &file_name, double time_seconds);

  /*!
   * Write a 1D array `input` to a variable `variable_name` in the file `file_name`.
   *
   * Data in `input` are written without modification.
   *
   * The array `input` is stored *redundantly* on all MPI ranks.
   *
   * FIXME: writing to the time variable will change the length of the time dimension.
   *
   * This version is used in Python scripts that create input files for testing. It should
   * not be used in PISM itself.
   */
  void write_array(const std::string &file_name, const std::string &variable_name,
                   const std::vector<unsigned int> &start, const std::vector<unsigned int> &count,
                   const std::vector<double> &input);

  /*!
   * Write a 1D array `input` to a variable `variable_name` in the file `file_name`,
   * converting from internal to "output" units if necessary.
   *
   * The array `input` is stored *redundantly* on all MPI ranks.
   *
   * FIXME: writing to the time variable will change the length of the time dimension.
   */
  void write_array(const std::string &file_name, const VariableMetadata &metadata,
                   const std::vector<unsigned int> &start, const std::vector<unsigned int> &count,
                   const std::vector<double> &input);

  /*!
   * Write a text (string) variable.
   *
   * The array `input` is stored *redundantly* on all MPI ranks.
   *
   */
  void write_text(const std::string &file_name, const std::string &variable_name,
                  const std::vector<unsigned int> &start, const std::vector<unsigned int> &count,
                  const std::string &input);

  /*!
   * Write a 2D or 3D array `input` described by `metadata` to the file `file_name`.
   *
   * Write coordinate variables (`x`, `y`, `z`, etc) required by this variable.
   *
   * Convert from internal to output units, if necessary.
   *
   * May be a no-op if this variable is time-independent and was written already.
   *
   * The array `input` is distributed across MPI ranks in the communicator used to create
   * this `OutputWriter` instance. Uses domain decomposition information provided to
   * `define_spatial_variable()`.
   */
  void write_spatial_variable(const std::string &file_name, const std::string &variable_name,
                              const double *input);

  /*!
   * Write a scalar time-dependent variable.
   */
  void write_timeseries_variable(const std::string &file_name, const VariableMetadata &metadata,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<double> &input);
  /*!
   * Indicate that the file `file_name` should be open for appending.
   *
   * This implies that if `file_name` should not be deleted if it already exists.
   *
   * May require reading time dimension length and the last value of time from `file_name`.
   */
  void append(const std::string &file_name);

  /*!
   * Ensure that all requested write operations are complete.
   *
   * May be a no-op in the context of asynchronous writing.
   */
  void sync(const std::string &file_name);

  /*!
   * Possibly close the file `file_name`. Used to indicate that the user does not intend
   * to write to the file `file_name` any more.
   */
  void close(const std::string &file_name);

  /*!
   * Return the length of the time dimension (possibly cached to avoid reading from the
   * file or communication).
   */
  unsigned int time_dimension_length(const std::string &file_name);

  /*!
   * Return the last value of the coordinate variable "time" (possibly cached to avoid
   * reading from the file or communication).
   *
   * Used when appending to an existing file.
   */
  double last_time_value(const std::string &file_name);
  
protected:
  /*!
   * Return the MPI communicator
   */
  MPI_Comm comm() const;

  /*!
   * Return the domain decomposition information for the variable `variable_name`.
   */
  const grid::DistributedGridInfo &grid_info(const std::string &variable_name) const;

  /*!
   * Return the metadata for the spatial variable `variable_name`.
   */
  const SpatialVariableMetadata &spatial_variable_info(const std::string &variable_name) const;

  /*!
   * Return `true` if variable `variable_name` was already written to the file
   * `file_name`. Used to avoid writing coordinate variables and time-independent 2D and
   * 3D arrays more than once.
   */
  bool &already_written(const std::string &file_name, const std::string &variable_name,
                        bool time_dependent);

  /*!
   * Return the name of the time dimension and the corresponding coordinate variable.
   */
  const std::string &time_name() const;

  /*!
   * Implementation of set_global_attributes()
   */
  virtual void
  set_global_attributes_impl(const std::string &file_name,
                             const std::map<std::string, std::string> &strings,
                             const std::map<std::string, std::vector<double> > &numbers) = 0;

  /*!
   * Implementation of define_dimension()
   */
  virtual void define_dimension_impl(const std::string &file_name, const std::string &name,
                                     unsigned int length) = 0;

  /*!
   * Implementation of define_variable()
   */
  virtual void define_variable_impl(const std::string &file_name, const VariableMetadata &metadata,
                                    const std::vector<std::string> &dims) = 0;

  /*!
   * Implementation of append_time()
   */
  virtual void append_time_impl(const std::string &file_name, double time_seconds) = 0;

  /*!
   * Implementation of append_history()
   */
  virtual void append_history_impl(const std::string &file_name, const std::string &text) = 0;

  /*!
   * Implementation of time_dimension_length()
   */
  virtual unsigned int time_dimension_length_impl(const std::string &file_name) = 0;

  /*!
   * Implementation of last_time_value()
   */
  virtual double last_time_value_impl(const std::string &file_name) = 0;

  /*!
   * Implementation of write_array()
   */
  virtual void write_array_impl(const std::string &file_name, const std::string &variable_name,
                                const std::vector<unsigned int> &start,
                                const std::vector<unsigned int> &count, const double *data) = 0;

  /*!
   * Implementation of write_text()
   */
  virtual void write_text_impl(const std::string &file_name, const std::string &variable_name,
                               const std::vector<unsigned int> &start,
                               const std::vector<unsigned int> &count,
                               const std::string &input) = 0;

  /*!
   * Implementation of write_spatial_variable()
   */
  virtual void write_spatial_variable_impl(const std::string &file_name,
                                           const std::string &variable_name,
                                           const double *data) = 0;

  /*!
   * Implementation of append()
   */
  virtual void append_impl(const std::string &file_name) = 0;

  /*!
   * Implementation of sync()
   */
  virtual void sync_impl(const std::string &file_name) = 0;

  /*!
   * Implementation of close()
   */
  virtual void close_impl(const std::string &file_name) = 0;

  const std::string &experiment_id() const;
private:
  void define_experiment_id(const std::string &file_name,
                            std::shared_ptr<units::System> unit_system);
  void write_experiment_id(const std::string &file_name);

  struct Impl;
  Impl *m_impl;
};

/*!
 * Wrapper class used to make OutputWriter a bit easier to use.
 *
 * See documentation of OutputWriter for details.
 *
 * Does not open the file when created, allowing one to call `append()` to indicate that a
 * file should not be over-written.
 */
class OutputFile {
public:
  OutputFile(std::shared_ptr<OutputWriter> writer, const std::string &file_name);

  void add_extra_attributes(const std::map<std::string, std::string> &attributes) const;

  void define_dimension(const std::string &dimension_name, unsigned int length) const;

  void define_variable(const VariableMetadata &metadata, const std::vector<std::string> &dims) const;

  void define_spatial_variable(const SpatialVariableMetadata &metadata,
                               const grid::DistributedGridInfo &grid) const;

  void define_timeseries_variable(const VariableMetadata &metadata) const;

  void set_global_attributes(const std::map<std::string, std::string> &strings,
                             const std::map<std::string, std::vector<double> > &numbers) const;

  void append_time(double time_seconds) const;

  void append_history(const std::string &text) const;

  void write_array(const std::string &variable_name, const std::vector<unsigned int> &start,
                   const std::vector<unsigned int> &count, const std::vector<double> &input) const;

  void write_array(const VariableMetadata &metadata, const std::vector<unsigned int> &start,
                   const std::vector<unsigned int> &count, const std::vector<double> &input) const;

  void write_spatial_variable(const std::string &variable_name, const double *input) const;

  void write_timeseries_variable(const VariableMetadata &metadata,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<double> &input) const;

  void write_text(const std::string &variable_name, const std::vector<unsigned int> &start,
                  const std::vector<unsigned int> &count, const std::string &input) const;

  void append();

  void sync();

  void close();

  unsigned int time_dimension_length() const;

  double last_time_value() const;

  std::string name() const;
private:
  std::string m_file_name;
  std::shared_ptr<OutputWriter> m_writer;
};

} // namespace pism

#endif /* PISM_OUTPUTWRITER_H */
