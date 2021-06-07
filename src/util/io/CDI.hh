// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019, 2021 PISM Authors
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

#include "NCFile.hh"
#include <map>
#include <set>

namespace pism {
namespace io {

class CDI : public NCFile {
public:
  CDI(MPI_Comm com);
  virtual ~CDI();

protected:
  void open_impl(const std::string &filename, IO_Mode mode,
                 //const std::map<std::string, int> &varsi = std::map<std::string, int>(),
                 int FileID = -1, const std::map<std::string, int> &dimsa = std::map<std::string, int>());

  void create_impl(const std::string &filename, int FileID = -1, int filetype = 0);

  void sync_impl() const;

  void close_impl();

  // redef/enddef
  void enddef_impl() const;

  void redef_impl() const;

  // dim
  void def_dim_impl(const std::string &name, size_t length, int dim) const;

  void inq_dimid_impl(const std::string &dimension_name, bool &exists) const;

  void inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const;

  void inq_unlimdim_impl(std::string &result) const;

  // var
  void def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const;

  void get_vara_double_impl(const std::string &variable_name, const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count, double *ip) const;

  void put_vara_double_impl(const std::string &variable_name, const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count, const double *op) const;

  void get_varm_double_impl(const std::string &variable_name, const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count, const std::vector<unsigned int> &imap,
                            double *ip) const;

  void write_darray_impl(const std::string &variable_name, const IceGrid &grid, unsigned int z_count,
                         unsigned int record, const double *input);

  void inq_nvars_impl(int &result) const;

  void inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const;

  void inq_varnatts_impl(const std::string &variable_name, int &result) const;

  void inq_varid_impl(const std::string &variable_name, bool &exists) const;

  void inq_varname_impl(unsigned int j, std::string &result) const;

  // att
  void get_att_double_impl(const std::string &variable_name, const std::string &att_name,
                           std::vector<double> &result) const;

  void get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const;

  void put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type xtype,
                           const std::vector<double> &data) const;

  void put_att_text_impl(const std::string &variable_name, const std::string &att_name, const std::string &value) const;

  void inq_attname_impl(const std::string &variable_name, unsigned int n, std::string &result) const;

  void inq_atttype_impl(const std::string &variable_name, const std::string &att_name, IO_Type &result) const;

  // misc
  void set_fill_impl(int fillmode, int &old_modep) const;

  void del_att_impl(const std::string &variable_name, const std::string &att_name) const;

  // new functions
  void create_grid_impl(int lengthx, int lengthy) const;
  void define_timestep_impl(int tsID) const;
  void def_ref_date_impl(double time) const;
  std::map<std::string, int> get_var_map_impl();
  std::map<std::string, int> get_dim_map_impl();
  void def_vlist_impl() const;
  void set_diagvars_impl(const std::set<std::string> &variables) const;
  void set_bdiag_impl(bool value) const;
  int get_ncstreamID_impl() const;
  int get_ncvlistID_impl() const;
  void set_calendar_impl(double year_length, const std::string &calendar_string) const;

private:
  // main
  mutable int m_gridID;
  mutable int m_tID;
  mutable int m_vlistID;
  mutable std::map<std::string, int> m_varsID;
  mutable std::map<std::string, int> m_dimsAxis;
  mutable std::map<std::string, int> m_zID;
  mutable std::set<std::string> m_diagvars;

  // auxiliary
  mutable int m_gridsID;
  mutable int m_zsID;
  mutable int m_gridgID;

  //switch
  mutable bool m_beforediag;

  // calendar
  mutable double m_year_length;
  mutable std::string m_calendar_string;
  mutable int m_cdi_calendar;
  mutable double m_days_year;

  // initialize class for opened file
  void map_varsID() const;
  void map_zaxisID() const;

  // file definition helpers
  void def_vlist() const;
  void def_zs() const;

  // inquire current timestep helper
  int inq_current_timestep() const;

  // calendar conversion
  double year_calendar(double time) const;
  long int day_calendar(double nyearsf) const;
  void monthday_calendar(int year, int doy, int *month, int *day) const;

  // inquire attribute helper
  void inq_att_impl(int varID, int attnum, char *attname, int *atttype, int *attlen) const;

  // inquire dimension length wrappers
  unsigned int inq_dimlen_t(const std::string &dimension_name) const;
  unsigned int inq_dimlen_x(const std::string &dimension_name) const;
  unsigned int inq_dimlen_y(const std::string &dimension_name) const;
  unsigned int inq_dimlen_z(const std::string &dimension_name) const;
  void wrapup_inq_dimlen() const;
  typedef unsigned int (CDI::*pInqDimlen)(const std::string &) const;
  mutable std::vector<pInqDimlen> pvcInqDimlen;

  // put attribute text wrappers
  void put_att_text_units_x_impl(const std::string &variable_name, const std::string &value) const;
  void put_att_text_longname_x_impl(const std::string &variable_name, const std::string &value) const;
  void put_att_text_units_y_impl(const std::string &variable_name, const std::string &value) const;
  void put_att_text_longname_y_impl(const std::string &variable_name, const std::string &value) const;
  void put_att_text_units_z_impl(const std::string &variable_name, const std::string &value) const;
  void put_att_text_longname_z_impl(const std::string &variable_name, const std::string &value) const;
  void wrapup_put_att_text() const;
  typedef void (CDI::*pPutAttT)(const std::string &, const std::string &) const;
  mutable std::vector<std::vector<pPutAttT> > pvcPutAttT;
  mutable std::map<std::string, int> m_DimAtt;

  // define variable wrappers
  void def_var_scalar_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const;
  void def_var_mscalar_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const;
  void def_var_multi_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const;
};
} // namespace io
} // namespace pism
