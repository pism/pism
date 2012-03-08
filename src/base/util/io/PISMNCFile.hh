// Copyright (C) 2012 PISM Authors
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

#ifndef _PISMNCWRAPPER_H_
#define _PISMNCWRAPPER_H_

#include <mpi.h>
#include <string>
#include <vector>

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>		// nc_type
// Note: as far as I (CK) can tell, MPI_INCLUDED is a MPICH invention.

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

//! \brief The PISM wrapper for a subset of the NetCDF C API.
/*!
 * The goal of this class is to hide the fact that we need to communicate data
 * to and from the processor zero. Using this wrapper we should be able to
 * write code that looks good and works both on 1-processor and
 * multi-processor systems.
 *
 * Moreover, this way we can switch underlying I/O implementations.
 *
 * Notes:
 * - It uses C++ STL strings instead of C character arrays
 * - It hides NetCDF ncid, dimid and varid and uses strings to reference
 *   dimensions and variables instead.
 * - This class does not and should not use any PETSc API calls.
 * - This wrapper provides access to a very small portion of the NetCDF C API.
 *   (Only calls used in PISM.) This is intentional.
 * - Methods of this class should do what corresponding NetCDF C API calls do,
 *   no more and no less.
 */
class PISMNCFile
{
public:
  PISMNCFile(MPI_Comm com, int rank);
  virtual ~PISMNCFile();

  // open/create/close
  virtual int open(string filename, int mode) = 0;

  virtual int create(string filename) = 0;

  virtual int close() = 0;

  // redef/enddef
  virtual int enddef() const = 0;

  virtual int redef() const = 0;

  // dim
  virtual int def_dim(string name, size_t length) const = 0;

  virtual int inq_dimid(string dimension_name, bool &exists) const = 0;

  virtual int inq_dimlen(string dimension_name, unsigned int &result) const = 0;

  virtual int inq_unlimdim(string &result) const = 0;

  // var
  virtual int def_var(string name, nc_type nctype, vector<string> dims) const = 0;

  virtual int get_vara_double(string variable_name,
                              vector<unsigned int> start,
                              vector<unsigned int> count,
                              double *ip) const = 0;

  virtual int put_vara_double(string variable_name,
                              vector<unsigned int> start,
                              vector<unsigned int> count,
                              const double *op) const = 0;

  virtual int get_varm_double(string variable_name,
                              vector<unsigned int> start,
                              vector<unsigned int> count,
                              vector<unsigned int> imap, double *ip) const = 0;

  virtual int put_varm_double(string variable_name,
                              vector<unsigned int> start,
                              vector<unsigned int> count,
                              vector<unsigned int> imap, const double *op) const = 0;

  virtual int inq_nvars(int &result) const = 0;

  virtual int inq_vardimid(string variable_name, vector<string> &result) const = 0;

  virtual int inq_varnatts(string variable_name, int &result) const = 0;

  virtual int inq_varid(string variable_name, bool &exists) const = 0;

  virtual int inq_varname(unsigned int j, string &result) const = 0;

  // att
  virtual int get_att_double(string variable_name, string att_name, vector<double> &result) const = 0;

  virtual int get_att_text(string variable_name, string att_name, string &result) const = 0;

  virtual int put_att_double(string variable_name, string att_name, nc_type xtype, vector<double> &data) const = 0;

  virtual int put_att_double(string variable_name, string att_name, nc_type xtype, double value) const;

  virtual int put_att_text(string variable_name, string att_name, string value) const = 0;

  virtual int inq_attname(string variable_name, unsigned int n, string &result) const = 0;

  virtual int inq_atttype(string variable_name, string att_name, nc_type &result) const = 0;

  // misc
  virtual int set_fill(int fillmode, int &old_modep) const = 0;

  string get_filename() const;

protected:

  void check(int return_code) const;

  int rank;
  MPI_Comm com;

  int ncid;
  string filename;
  mutable bool define_mode;
};

#endif /* _PISMNCWRAPPER_H_ */
