// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019 PISM Authors
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

class CDI : public NCFile
{
public:
	CDI(MPI_Comm com);
	virtual ~CDI();
protected:
	void open_impl(const std::string &filename, IO_Mode mode, const std::map<std::string, int> &varsi = std::map<std::string, int>(),
                       int FileID = -1);

	void create_impl(const std::string &filename, int FileID = -1);

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
	void def_var_impl(const std::string &name, IO_Type nctype,
              const std::vector<std::string> &dims) const;

	void get_vara_double_impl(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      double *ip) const;

	void put_vara_double_impl(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      const double *op) const;

	void get_varm_double_impl(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      const std::vector<unsigned int> &imap,
                      double *ip) const;

	void write_darray_impl(const std::string &variable_name,
                      const IceGrid &grid,
                      unsigned int z_count,
                      unsigned int record,
                      const double *input);

	void inq_nvars_impl(int &result) const;

	void inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const;

	void inq_varnatts_impl(const std::string &variable_name, int &result) const;

	void inq_varid_impl(const std::string &variable_name, bool &exists) const;

	void inq_varname_impl(unsigned int j, std::string &result) const;

	// att
	void get_att_double_impl(const std::string &variable_name, const std::string &att_name,
                     std::vector<double> &result) const;

	void get_att_text_impl(const std::string &variable_name, const std::string &att_name,
                   std::string &result) const;

	void put_att_double_impl(const std::string &variable_name, const std::string &att_name,
                     IO_Type xtype, const std::vector<double> &data) const;

	void put_att_text_impl(const std::string &variable_name, const std::string &att_name,
                   const std::string &value) const;

	void inq_attname_impl(const std::string &variable_name, unsigned int n,
                  std::string &result) const;

	void inq_atttype_impl(const std::string &variable_name, const std::string &att_name,
                  IO_Type &result) const;

	// misc
	void set_fill_impl(int fillmode, int &old_modep) const;

	void del_att_impl(const std::string &variable_name, const std::string &att_name) const;

	// new functions
	void create_grid_impl(int lengthx, int lengthy) const;
	void define_timestep_impl(int tsID) const;
	void def_ref_date_impl(double time) const;
	std::map<std::string, int> get_var_map_impl();
        void def_vlist_impl() const;
	void set_diagvars_impl(const std::set<std::string> &variables) const;
	void set_bdiag_impl(bool value) const;
	void set_ncgridIDs_impl(const std::vector<int>& gridIDs) const;
	std::vector<int> get_ncgridIDs_impl() const;
	int get_ncstreamID_impl() const;
	int get_ncvlistID_impl() const;

private:
	mutable int m_gridID;
	mutable int m_gridsID;
	mutable int m_zID;
	mutable int m_zbID;
	mutable int m_zsID;
	mutable int m_tID;
	mutable int m_vlistID;
	mutable std::map<std::string,int> m_varsID;
    mutable std::vector<std::string> m_dims_name;
	mutable bool m_firststep;
	mutable bool m_beforediag;
	mutable std::set<std::string> m_diagvars;
	mutable bool m_gridexist;
	mutable bool m_istimedef;
	void def_var_scalar_impl(const std::string &name, IO_Type nctype,
	      const std::vector<std::string> &dims) const;
    void def_var_multi_impl(const std::string &name, IO_Type nctype,
              const std::vector<std::string> &dims) const;
    int inq_current_timestep() const;
    void put_att_text_dims_impl(const std::string &variable_name, const std::string &att_name,
                   const std::string &value) const;
	void put_att_text_x_impl(const std::string &att_name, const std::string &value) const;
	void put_att_text_y_impl(const std::string &att_name, const std::string &value) const;
	void put_att_text_z_impl(int zaxisID, const std::string &att_name,
                                const std::string &value) const;
	void def_main_dim(const std::string &name, size_t length) const;
	double year_gregorian(double time) const;
	long int day_gregorian(double nyearsf) const;
	void put_dims_double(const std::string &variable_name, const double *op) const;
	void map_varsID() const;
	void inq_att_impl(int varID, int attnum, char* attname, int *atttype, int *attlen) const;
	mutable std::map<std::string,int> m_dimsAxis;
	mutable std::map<std::string,int> m_zID;

	// functions wrappers
	void def_x_dim(const std::string &name, size_t length) const;
	void def_y_dim(const std::string &name, size_t length) const;
	void def_z_dim(const std::string &name, size_t length) const;
	void def_t_dim(const std::string &name, size_t length) const;
	void wrapup_def_dim() const;
	typedef void(CDI::*pDefDim)(const std::string&, size_t) const;
	mutable std::vector<pDefDim> pvcDefDim;
	unsigned int CDI::inq_dimlen_t(const std::string &dimension_name) const;
	unsigned int CDI::inq_dimlen_x(const std::string &dimension_name) const;
	unsigned int CDI::inq_dimlen_y(const std::string &dimension_name) const;
	unsigned int CDI::inq_dimlen_z(const std::string &dimension_name) const;
	void CDI::wrapup_inq_dimlen() const;
	typedef unsigned int(CDI::*pInqDimlen)(const std::string&) const;
	mutable std::vector<pInqDimlen> pvcInqDimlen;
	
};
}
}
