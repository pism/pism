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

#ifndef _PISMNC3FILE_H_
#define _PISMNC3FILE_H_

#include "PISMNCFile.hh"

class PISMNC3File : public PISMNCFile
{
public:
  PISMNC3File(MPI_Comm com, int rank);
  virtual ~PISMNC3File();

  // open/create/close
  int open(string filename, int mode);

  int create(string filename);

  int close();

  // redef/enddef
  int enddef() const;

  int redef() const;

  // dim
  int def_dim(string name, size_t length) const;

  int inq_dimid(string dimension_name, bool &exists) const;

  int inq_dimlen(string dimension_name, unsigned int &result) const;

  int inq_unlimdim(string &result) const;

  int inq_dimname(int j, string &result) const;

  int inq_ndims(int &result) const;

  // var
  int def_var(string name, PISM_IO_Type nctype, vector<string> dims) const;

  int get_vara_double(string variable_name,
                      vector<unsigned int> start,
                      vector<unsigned int> count,
                      double *ip) const;

  int put_vara_double(string variable_name,
                      vector<unsigned int> start,
                      vector<unsigned int> count,
                      double *op) const;

  int get_varm_double(string variable_name,
                      vector<unsigned int> start,
                      vector<unsigned int> count,
                      vector<unsigned int> imap, double *ip) const;

  int put_varm_double(string variable_name,
                      vector<unsigned int> start,
                      vector<unsigned int> count,
                      vector<unsigned int> imap, double *op) const;

  int inq_nvars(int &result) const;

  int inq_vardimid(string variable_name, vector<string> &result) const;

  int inq_varnatts(string variable_name, int &result) const;

  int inq_varid(string variable_name, bool &exists) const;

  int inq_varname(unsigned int j, string &result) const;

  int inq_vartype(string variable_name, PISM_IO_Type &result) const;
  // att
  int get_att_double(string variable_name, string att_name, vector<double> &result) const;

  int get_att_text(string variable_name, string att_name, string &result) const;

  using PISMNCFile::put_att_double;
  int put_att_double(string variable_name, string att_name, PISM_IO_Type xtype, vector<double> &data) const;

  int put_att_text(string variable_name, string att_name, string value) const;

  int inq_attname(string variable_name, unsigned int n, string &result) const;

  int inq_atttype(string variable_name, string att_name, PISM_IO_Type &result) const;

  // misc
  int set_fill(int fillmode, int &old_modep) const;

  virtual string get_format() const;
private:
  int get_var_double(string variable_name,
                     vector<unsigned int> start,
                     vector<unsigned int> count,
                     vector<unsigned int> imap, double *ip,
                     bool mapped) const;

  int put_var_double(string variable_name,
                     vector<unsigned int> start,
                     vector<unsigned int> count,
                     vector<unsigned int> imap, const double *op,
                     bool mapped) const;
};

#endif /* _PISMNC3FILE_H_ */
