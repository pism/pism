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

#ifndef _PISMPNCFILE_H_
#define _PISMPNCFILE_H_

#include "PISMNCFile.hh"

namespace pism {

//! \brief PISM's PnetCDF I/O wrapper.
class PNCFile : public NCFile
{
public:
  PNCFile(MPI_Comm com);
  virtual ~PNCFile();

  // open/create/close
  int open(std::string filename, IO_Mode mode);

  int create(std::string filename);

  int close();

  // redef/enddef
  int enddef() const;

  int redef() const;

  // dim
  int def_dim(std::string name, size_t length) const;

  int inq_dimid(std::string dimension_name, bool &exists) const;

  int inq_dimlen(std::string dimension_name, unsigned int &result) const;

  int inq_unlimdim(std::string &result) const;

  int inq_dimname(int j, std::string &result) const;

  int inq_ndims(int &result) const;

  // var
  int def_var(std::string name, IO_Type nctype, std::vector<std::string> dims) const;

  int get_vara_double(std::string variable_name,
                      std::vector<unsigned int> start,
                      std::vector<unsigned int> count,
                      double *ip) const;

  int put_vara_double(std::string variable_name,
                      std::vector<unsigned int> start,
                      std::vector<unsigned int> count,
                      const double *op) const;

  int get_varm_double(std::string variable_name,
                      std::vector<unsigned int> start,
                      std::vector<unsigned int> count,
                      std::vector<unsigned int> imap,
                      double *ip) const;

  int put_varm_double(std::string variable_name,
                      std::vector<unsigned int> start,
                      std::vector<unsigned int> count,
                      std::vector<unsigned int> imap,
                      const double *op) const;

  int inq_nvars(int &result) const;

  int inq_vardimid(std::string variable_name, std::vector<std::string> &result) const;

  int inq_varnatts(std::string variable_name, int &result) const;

  int inq_varid(std::string variable_name, bool &exists) const;

  int inq_varname(unsigned int j, std::string &result) const;

  int inq_vartype(std::string variable_name, IO_Type &result) const;

  // att
  int get_att_double(std::string variable_name, std::string att_name, std::vector<double> &result) const;

  int get_att_text(std::string variable_name, std::string att_name, std::string &result) const;

  using NCFile::put_att_double;
  int put_att_double(std::string variable_name, std::string att_name, IO_Type xtype, const std::vector<double> &data) const;

  int put_att_text(std::string variable_name, std::string att_name, std::string value) const;

  int inq_attname(std::string variable_name, unsigned int n, std::string &result) const;

  int inq_atttype(std::string variable_name, std::string att_name, IO_Type &result) const;

  // misc
  int set_fill(int fillmode, int &old_modep) const;

  virtual std::string get_format() const;

  std::vector<std::string> mpi_io_hints;
protected:
  virtual int integer_open_mode(IO_Mode input) const;
  void check(int return_code) const;

private:
  int get_var_double(std::string variable_name,
                     std::vector<unsigned int> start,
                     std::vector<unsigned int> count,
                     std::vector<unsigned int> imap, double *ip,
                     bool mapped) const;

  int put_var_double(std::string variable_name,
                     std::vector<unsigned int> start,
                     std::vector<unsigned int> count,
                     std::vector<unsigned int> imap, const double *op,
                     bool mapped) const;

  void init_hints();

  MPI_Info mpi_info;            // MPI hints
};

} // end of namespace pism

#endif /* _PISMPNCFILE_H_ */
