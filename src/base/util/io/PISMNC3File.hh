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

#ifndef _PISMNC3FILE_H_
#define _PISMNC3FILE_H_

#include "PISMNCFile.hh"

namespace pism {

class NC3File : public NCFile
{
public:
  NC3File(MPI_Comm com);
  virtual ~NC3File();

protected:
  // implementations:
  // open/create/close
  int open_impl(const std::string &filename, IO_Mode mode);

  int create_impl(const std::string &filename);

  int close_impl();

  // redef/enddef
  int enddef_impl() const;

  int redef_impl() const;

  // dim
  int def_dim_impl(const std::string &name, size_t length) const;

  int inq_dimid_impl(const std::string &dimension_name, bool &exists) const;

  int inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const;

  int inq_unlimdim_impl(std::string &result) const;

  int inq_dimname_impl(int j, std::string &result) const;

  int inq_ndims_impl(int &result) const;

  // var
  int def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const;

  int get_vara_double_impl(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      double *ip) const;

  int put_vara_double_impl(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      const double *op) const;

  int get_varm_double_impl(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      const std::vector<unsigned int> &imap,
                      double *ip) const;

  int put_varm_double_impl(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      const std::vector<unsigned int> &imap,
                      const double *op) const;

  int inq_nvars_impl(int &result) const;

  int inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const;

  int inq_varnatts_impl(const std::string &variable_name, int &result) const;

  int inq_varid_impl(const std::string &variable_name, bool &exists) const;

  int inq_varname_impl(unsigned int j, std::string &result) const;

  int inq_vartype_impl(const std::string &variable_name, IO_Type &result) const;
  // att
  int get_att_double_impl(const std::string &variable_name, const std::string &att_name, std::vector<double> &result) const;

  int get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const;

  using NCFile::put_att_double;
  int put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type xtype, const std::vector<double> &data) const;

  int put_att_text_impl(const std::string &variable_name, const std::string &att_name, const std::string &value) const;

  int inq_attname_impl(const std::string &variable_name, unsigned int n, std::string &result) const;

  int inq_atttype_impl(const std::string &variable_name, const std::string &att_name, IO_Type &result) const;

  // misc
  int set_fill_impl(int fillmode, int &old_modep) const;

  std::string get_format_impl() const;
private:
  int m_rank;
  int integer_open_mode(IO_Mode input) const;

  int get_var_double(const std::string &variable_name,
                     const std::vector<unsigned int> &start,
                     const std::vector<unsigned int> &count,
                     const std::vector<unsigned int> &imap, double *ip,
                     bool mapped) const;

  int put_var_double(const std::string &variable_name,
                     const std::vector<unsigned int> &start,
                     const std::vector<unsigned int> &count,
                     const std::vector<unsigned int> &imap, const double *op,
                     bool mapped) const;
};

} // end of namespace pism

#endif /* _PISMNC3FILE_H_ */
