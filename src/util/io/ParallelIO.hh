/* Copyright (C) 2019, 2020 PISM Authors
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

#include "NCFile.hh"
#include "IO_Flags.hh"

namespace pism {
namespace io {

class ParallelIO : public NCFile {
public:
  ParallelIO(MPI_Comm com, int iosysid, IO_Backend iotype);
  virtual ~ParallelIO();

  static IO_Backend best_iotype(bool netcdf3);
protected:
  // open/create/close
  void open_impl(const std::string &filename, IO_Mode mode, const std::map<std::string, int> &varsi = std::map<std::string, int>(),
                 int FileID = -1);
  void create_impl(const std::string &filename, int FileID = -1);
  void sync_impl() const;
  void close_impl();

  void set_compression_level_impl(int level) const;

  // redef/enddef
  void enddef_impl() const;

  void redef_impl() const;

  // dim
  void def_dim_impl(const std::string &name, size_t length, int dim) const;

  void inq_dimid_impl(const std::string &dimension_name, bool &exists) const;

  void inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const;

  void inq_unlimdim_impl(std::string &result) const;

  // var
  void def_var_impl(const std::string &name, IO_Type nctype,
                   const std::vector<std::string> &dims) const;

  void def_var_chunking_impl(const std::string &name,
                            std::vector<size_t> &dimensions) const;

  void get_vara_double_impl(const std::string &variable_name,
                           const std::vector<unsigned int> &start,
                           const std::vector<unsigned int> &count,
                           double *ip) const;

  void put_vara_double_impl(const std::string &variable_name,
                           const std::vector<unsigned int> &start,
                           const std::vector<unsigned int> &count,
                           const double *op) const;

  void write_darray_impl(const std::string &variable_name,
                         const IceGrid &grid,
                         unsigned int z_count,
                         unsigned int record,
                         const double *input);

  void get_varm_double_impl(const std::string &variable_name,
                           const std::vector<unsigned int> &start,
                           const std::vector<unsigned int> &count,
                           const std::vector<unsigned int> &imap,
                           double *ip) const;

  void inq_nvars_impl(int &result) const;

  void inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const;

  void inq_varnatts_impl(const std::string &variable_name, int &result) const;

  void inq_varid_impl(const std::string &variable_name, bool &exists) const;

  void inq_varname_impl(unsigned int j, std::string &result) const;

  // att
  void get_att_double_impl(const std::string &variable_name, const std::string &att_name, std::vector<double> &result) const;

  void get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const;

  void put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type xtype, const std::vector<double> &data) const;

  void put_att_text_impl(const std::string &variable_name, const std::string &att_name, const std::string &value) const;

  void inq_attname_impl(const std::string &variable_name, unsigned int n, std::string &result) const;

  void inq_atttype_impl(const std::string &variable_name, const std::string &att_name, IO_Type &result) const;

  // misc
  void set_fill_impl(int fillmode, int &old_modep) const;

  void del_att_impl(const std::string &variable_name, const std::string &att_name) const;

  // new functions empty because of CDI class
  void create_grid_impl(int lengthx, int lengthy) const;
  void define_timestep_impl(int tsID) const;
  void def_ref_date_impl(double time) const;

private:
  int m_iosysid;
  int m_iotype;
  int get_varid(const std::string &variable_name) const;
};

} // end of namespace io
} // end of namespace pism
