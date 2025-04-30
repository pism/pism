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
 *
 *    
 */
class OutputWriter {
public:
  OutputWriter(MPI_Comm comm, const Config &config);
  virtual ~OutputWriter();

  void add_extra_attributes(const std::string &file_name,
                            const std::map<std::string, std::string> &attributes);

  void define_dimension(const std::string &file_name, const std::string &dimension_name,
                        unsigned int length);

  void define_variable(const std::string &file_name, const VariableMetadata &metadata,
                       const std::vector<std::string> &dims);

  void define_spatial_variable(const std::string &file_name,
                               const SpatialVariableMetadata &metadata,
                               const grid::DistributedGridInfo &grid);

  void write_attributes(const std::string &file_name, const VariableMetadata &variable);

  void append_history(const std::string &file_name, const std::string &text);

  void append_time(const std::string &file_name, double time_seconds);

  void write_array(const std::string &file_name, const std::string &variable_name,
                   const std::vector<unsigned int> &start, const std::vector<unsigned int> &count,
                   const std::vector<double> &input);

  void write_array(const std::string &file_name, const VariableMetadata &metadata,
                   const std::vector<unsigned int> &start, const std::vector<unsigned int> &count,
                   const std::vector<double> &input);

  void write_spatial_variable(const std::string &file_name, const SpatialVariableMetadata &metadata,
                              const double *input);

  void append(const std::string &file_name);

  void sync(const std::string &file_name);

  void close(const std::string &file_name);

  unsigned int time_dimension_length(const std::string &file_name);

  double last_time_value(const std::string &file_name);
  
protected:
  MPI_Comm comm() const;

  const grid::DistributedGridInfo &grid_info(const std::string &variable_name) const;

  bool &already_written(const std::string &file_name, const std::string &variable_name);

  const std::string &time_name() const;

  virtual void write_attributes(const std::string &file_name, const std::string &var_name,
                                const std::map<std::string, std::string> &strings,
                                const std::map<std::string, std::vector<double> > &numbers,
                                io::Type output_type) = 0;

  virtual void define_dimension_impl(const std::string &file_name, const std::string &name,
                                     unsigned int length) = 0;

  virtual void define_variable_impl(const std::string &file_name, const VariableMetadata &metadata,
                                    const std::vector<std::string> &dims) = 0;

  virtual void append_time_impl(const std::string &file_name, double time_seconds) = 0;

  virtual void append_history_impl(const std::string &file_name, const std::string &text) = 0;

  virtual unsigned int time_dimension_length_impl(const std::string &file_name) = 0;

  virtual double last_time_value_impl(const std::string &file_name) = 0;

  virtual void write_array_impl(const std::string &file_name, const std::string &variable_name,
                                const std::vector<unsigned int> &start,
                                const std::vector<unsigned int> &count, const double *data) = 0;

  virtual void write_distributed_array_impl(const std::string &file_name,
                                            const std::string &variable_name,
                                            const std::vector<unsigned int> &start,
                                            const std::vector<unsigned int> &count,
                                            const double *data) = 0;

  virtual void append_impl(const std::string &file_name) = 0;

  virtual void sync_impl(const std::string &file_name) = 0;

  virtual void close_impl(const std::string &file_name) = 0;

private:
  struct Impl;
  Impl *m_impl;
};

class OutputFile {
public:
  OutputFile(std::shared_ptr<OutputWriter> writer, const std::string &file_name);

  void add_extra_attributes(const std::map<std::string, std::string> &attributes) const;

  void define_dimension(const std::string &dimension_name, unsigned int length) const;

  void define_variable(const VariableMetadata &metadata, const std::vector<std::string> &dims) const;

  void define_spatial_variable(const SpatialVariableMetadata &metadata,
                               const grid::DistributedGridInfo &grid) const;

  void write_attributes(const VariableMetadata &variable) const;

  void append_time(double time_seconds) const;

  void append_history(const std::string &text) const;

  void write_array(const std::string &variable_name, const std::vector<unsigned int> &start,
                   const std::vector<unsigned int> &count, const std::vector<double> &input) const;

  void write_array(const VariableMetadata &metadata, const std::vector<unsigned int> &start,
                   const std::vector<unsigned int> &count, const std::vector<double> &input) const;

  void write_spatial_variable(const SpatialVariableMetadata &metadata, const double *input) const;

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
