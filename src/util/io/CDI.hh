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

namespace pism {
namespace io {

class CDI : public NCFile
{
public:
	CDI(MPI_Comm com);
	virtual ~CDI();
protected:
	void open_impl(const std::string &filename, IO_Mode mode);

	void create_impl(const std::string &filename);

	void sync_impl() const;

	void close_impl();

	// redef/enddef
	void enddef_impl() const;

	void redef_impl() const;

	// dim
	void def_dim_impl(const std::string &name, size_t length) const;

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

private:
	mutable int m_gridID;
	mutable int m_zID;
	mutable int m_zbID;
        mutable int m_zsID;
	mutable int m_tID;
	mutable int m_vlistID;
	mutable std::map<std::string, int> m_varsID;
};
}
}
