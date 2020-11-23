// Copyright (C) 2012, 2013, 2014, 2015, 2017, 2019, 2020 PISM Authors
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

#ifndef _PISMNC4FILE_H_
#define _PISMNC4FILE_H_

#include "NCFile.hh"

namespace pism {
namespace io {

class NC4File : public NCFile
{
public:
  NC4File(MPI_Comm com, unsigned int compression_level);
  virtual ~NC4File();

protected:
  // implementations:
  // open/create/close

  virtual void sync_impl() const;

  virtual void close_impl();

  // redef/enddef
  virtual void enddef_impl() const;

  virtual void redef_impl() const;

  // dim
  virtual void def_dim_impl(const std::string &name, size_t length) const;

  virtual void inq_dimid_impl(const std::string &dimension_name, bool &exists) const;

  virtual void inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const;

  virtual void inq_unlimdim_impl(std::string &result) const;

  // var
  virtual void def_var_chunking_impl(const std::string &name,
                                    std::vector<size_t> &dimensions) const;

  virtual void def_var_impl(const std::string &name,
                           IO_Type nctype, const std::vector<std::string> &dims) const;

  virtual void get_vara_double_impl(const std::string &variable_name,
                                   const std::vector<unsigned int> &start,
                                   const std::vector<unsigned int> &count,
                                   double *ip) const;

  virtual void put_vara_double_impl(const std::string &variable_name,
                                   const std::vector<unsigned int> &start,
                                   const std::vector<unsigned int> &count,
                                   const double *op) const;

  virtual void get_varm_double_impl(const std::string &variable_name,
                                   const std::vector<unsigned int> &start,
                                   const std::vector<unsigned int> &count,
                                   const std::vector<unsigned int> &imap, double *ip) const;

  virtual void inq_nvars_impl(int &result) const;

  virtual void inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const;

  virtual void inq_varnatts_impl(const std::string &variable_name, int &result) const;

  virtual void inq_varid_impl(const std::string &variable_name, bool &exists) const;

  virtual void inq_varname_impl(unsigned int j, std::string &result) const;

  // att
  virtual void get_att_double_impl(const std::string &variable_name, const std::string &att_name, std::vector<double> &result) const;

  virtual void get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const;

  virtual void put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type xtype, const std::vector<double> &data) const;

  virtual void put_att_text_impl(const std::string &variable_name, const std::string &att_name, const std::string &value) const;

  virtual void inq_attname_impl(const std::string &variable_name, unsigned int n, std::string &result) const;

  virtual void inq_atttype_impl(const std::string &variable_name, const std::string &att_name, IO_Type &result) const;

  // misc
  virtual void set_fill_impl(int fillmode, int &old_modep) const;

  virtual void del_att_impl(const std::string &variable_name, const std::string &att_name) const;

  // new functions empty because of CDI class
  virtual void create_grid_impl(int lengthx, int lengthy) const;
  virtual void define_timestep_impl(int tsID) const;
  virtual void write_timestep_impl() const;

protected:
  virtual void set_access_mode(int varid, bool mapped) const;
  virtual void get_put_var_double(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap, double *ip,
                                 bool get,
                                 bool mapped) const;

  mutable unsigned int m_compression_level;

  int get_varid(const std::string &variable_name) const;
};

} // end of namespace io
} // end of namespace pism

#endif /* _PISMNC4FILE_H_ */
