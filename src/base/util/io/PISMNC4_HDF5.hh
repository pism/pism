// Copyright (C) 2012, 2013, 2014 PISM Authors
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

#ifndef _NC4_HDF5_H_
#define _NC4_HDF5_H_

#include <hdf5.h>

#include "PISMNCFile.hh"

namespace pism {

class NC4_HDF5 : public NCFile {
public:
  NC4_HDF5(MPI_Comm com);
  virtual ~NC4_HDF5();

  // open/create/close
  virtual int open(std::string filename, IO_Mode mode);

  virtual int create(std::string filename);

  virtual int close();

  // redef/enddef
  virtual int enddef() const;

  virtual int redef() const;

  // dim
  virtual int def_dim(std::string name, size_t length) const;

  virtual int inq_dimid(std::string dimension_name, bool &exists) const;

  virtual int inq_dimlen(std::string dimension_name, unsigned int &result) const;

  virtual int inq_unlimdim(std::string &result) const;

  virtual int inq_dimname(int j, std::string &result) const;

  virtual int inq_ndims(int &result) const;

  // var
  virtual int def_var(std::string name, IO_Type nctype, std::vector<std::string> dims) const;

  virtual int get_vara_double(std::string variable_name,
                              std::vector<unsigned int> start,
                              std::vector<unsigned int> count,
                              double *ip) const;

  virtual int put_vara_double(std::string variable_name,
                              std::vector<unsigned int> start,
                              std::vector<unsigned int> count,
                              const double *op) const;

  virtual int get_varm_double(std::string variable_name,
                              std::vector<unsigned int> start,
                              std::vector<unsigned int> count,
                              std::vector<unsigned int> imap, double *ip) const;

  virtual int put_varm_double(std::string variable_name,
                              std::vector<unsigned int> start,
                              std::vector<unsigned int> count,
                              std::vector<unsigned int> imap, const double *op) const;

  virtual int inq_nvars(int &result) const;

  virtual int inq_vardimid(std::string variable_name, std::vector<std::string> &result) const;

  virtual int inq_varnatts(std::string variable_name, int &result) const;

  virtual int inq_varid(std::string variable_name, bool &exists) const;

  virtual int inq_varname(unsigned int j, std::string &result) const;

  virtual int inq_vartype(std::string variable_name, IO_Type &result) const;

  // att
  virtual int get_att_double(std::string variable_name, std::string att_name, std::vector<double> &result) const;

  virtual int get_att_text(std::string variable_name, std::string att_name, std::string &result) const;

  using NCFile::put_att_double;
  virtual int put_att_double(std::string variable_name, std::string att_name, IO_Type xtype,
                             const std::vector<double> &data) const;

  virtual int put_att_text(std::string variable_name, std::string att_name, std::string value) const;

  virtual int inq_attname(std::string variable_name, unsigned int n, std::string &result) const;

  virtual int inq_atttype(std::string variable_name, std::string att_name, IO_Type &result) const;

  // misc
  virtual int set_fill(int fillmode, int &old_modep) const;

  virtual std::string get_format() const;
protected:
  virtual int integer_open_mode(IO_Mode input) const;
  virtual void check(int return_code) const;

  hid_t file_id;
};

} // end of namespace pism

#endif /* _NC4_HDF5_H_ */
