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
#include <mpi.h>
#include <sstream>
#include <string.h>
#include <memory>
#include <algorithm>
#include <map>
 #include <iostream>

#include "CDI.hh"

#include "pism/util/io/pism_cdi_type_conversion.hh"
#include "pism/util/IceGrid.hh"

#include "pism/util/error_handling.hh"

extern "C"{
#include "cdipio.h"
#include "cdi.h"
#include "yaxt.h"
}

namespace pism {
namespace io {

CDI::CDI(MPI_Comm c) : NCFile(c) {
	m_vlistID = -1;
	m_varsID.clear();
	m_beforediag = true;
	m_gridexist = false;
	m_istimedef = false;
	wrapup_def_dim();
}

CDI::~CDI() {
}

void CDI::open_impl(const std::string &fname, IO_Mode mode, 
	                const std::map<std::string, int> &varsi, 
	                int FileID//,
	                //const std::map<std::string, int> &dimsa
	                ) {
	// the file is already created and opened - restore file info into the class
	m_file_id = FileID;
   	m_vlistID = streamInqVlist(m_file_id);
   	m_varsID = varsi;
   	//m_dimsAxis = dimsa;
	m_firststep = false;
}

void CDI::map_varsID() const {
	int Nvars, varID;
	char VarName[CDI_MAX_NAME];
	std::string VarString;
	Nvars = vlistNvars(m_vlistID);
	for (varID=0; varID<Nvars; varID++) {
		vlistInqVarName(m_vlistID, varID, VarName);
		VarString = VarName;
		m_varsID[VarString] = varID;
	}
}

void CDI::set_ncgridIDs_impl(const std::vector<int>& gridIDs) const {
	if (gridIDs.empty()) {
		m_gridID = -1;
		m_gridsID = -1;
		m_tID = -1;
                // m_zID scalar deleted
		m_zbID = -1;
		m_zsID = -1;		
	} else {
		m_gridID = gridIDs[0];
		m_gridsID = gridIDs[1];
		m_tID = gridIDs[2];
                // m_zID scalar deleted
		m_zbID = gridIDs[4];
		m_zsID = gridIDs[5];
		m_gridexist = true;
	}
}

std::vector<int> CDI::get_ncgridIDs_impl() const {
        // m_zID scalar deleted
	std::vector<int> gridIDs{m_gridID, m_gridsID, m_tID, m_zbID, m_zsID};
	return gridIDs;
}


void CDI::create_impl(const std::string &filename, int FileID) {
        if (FileID == -1) {
			m_file_id = streamOpenWrite(filename.c_str(), CDI_FILETYPE_NC2);
        	m_firststep = true;
        } else {
        	m_file_id = FileID;
        	m_firststep = false;
        }
}

void CDI::close_impl() {
	m_file_id = -1;
	m_varsID.clear();
	m_dimsAxis.clear();
	m_zID.clear();
}

// wrapup def_dim_impl function
void CDI::def_x_dim(const std::string &name, size_t length) const {
    gridDefXsize(m_gridID, (int)length);
    gridDefXname(m_gridID, name.c_str());
}

void CDI::def_y_dim(const std::string &name, size_t length) const {
    gridDefYsize(m_gridID, (int)length);
    gridDefYname(m_gridID, name.c_str());
}

void CDI::def_z_dim(const std::string &name, size_t length) const {
	// define z axis only if it's new
	if (!m_zID.count(name)) {
	    m_zID[name] = zaxisCreate(ZAXIS_GENERIC, (int)length);
    	zaxisDefName(m_zID[name],name.c_str());
    }
}

void CDI::def_t_dim(const std::string &name, size_t length) const {
	// define time axis if it was not done before
	if (m_tID == -1) {
		m_tID = taxisCreate(TAXIS_ABSOLUTE);
		taxisDefCalendar(m_tID, CALENDAR_365DAYS);
    	vlistDefTaxis(m_vlistID, m_tID);
	}
}

void CDI::wrapup_def_dim() const {
    pDefDim pFn = &CDI::def_x_dim;
	pvcDefDim.push_back(pFn);
	pFn = &CDI::def_y_dim;
	pvcDefDim.push_back(pFn);
	pFn = &CDI::def_z_dim;
	pvcDefDim.push_back(pFn);
	pFn = &CDI::def_t_dim;
	pvcDefDim.push_back(pFn);
}

void CDI::def_dim_impl(const std::string &name, size_t length, int dim) const {
	if (m_vlistID == -1) m_vlistID = vlistCreate();

	if (!m_gridexist) {
		if (dim != -1) {
			(this->*(pvcDefDim[dim]))(name, length);
			m_dimsAxis[name] = dim;
		}
		if (m_zsID == -1) m_zsID = zaxisCreate(ZAXIS_SURFACE, 1);
	} else {
		if (!m_istimedef) {
			if (m_firststep) vlistDefTaxis(m_vlistID, m_tID);
			m_istimedef = true;
		}
	}
}

// define reference date
// TODO support other calendars
double CDI::year_gregorian(double time) const {
	return -time / 365 / 24 / 60 / 60;
}

long int CDI::day_gregorian(double nyearsf) const {
	long int seconds = (nyearsf - (long int)nyearsf) * 86400;
        long int minutes, hours;
        hours = seconds / 3600;
        minutes = (seconds - (3600*hours)) / 60;
        seconds = (seconds - (3600*hours) - minutes*60);
	return hours * 10000 + minutes * 100 + seconds;
}

void CDI::def_ref_date_impl(double time) const {
	// Taking into account only gregorian calendar right now
	double nyearsf = year_gregorian(time);
	taxisDefVdate(m_tID, -(long int)nyearsf);
	long int daytime = day_gregorian(nyearsf);
	taxisDefVtime(m_tID, daytime);
}

// inquire if a dimension exists
void CDI::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
	exists = m_dimsAxis.count(dimension_name) ? true : false;
}

// wrapup inq_dimlen_impl function
unsigned int CDI::inq_dimlen_x(const std::string &dimension_name) const {
	return gridInqXsize(m_gridID);
}

unsigned int CDI::inq_dimlen_y(const std::string &dimension_name) const {
	return gridInqYsize(m_gridID);
}

unsigned int CDI::inq_dimlen_z(const std::string &dimension_name) const {
	return zaxisInqSize(m_zID[dimension_name]);
}

unsigned int CDI::inq_dimlen_t(const std::string &dimension_name) const {
	return inq_current_timestep() + 1;
}

void CDI::wrapup_inq_dimlen() const {
    pInqDimlen pFn = &CDI::inq_dimlen_x;
	pvcInqDimlen.push_back(pFn);
	pFn = &CDI::inq_dimlen_y;
	pvcInqDimlen.push_back(pFn);
	pFn = &CDI::inq_dimlen_z;
	pvcInqDimlen.push_back(pFn);
	pFn = &CDI::inq_dimlen_t;
	pvcInqDimlen.push_back(pFn);
}

void CDI::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
	int dim = m_dimsAxis[dimension_name];
	result = (this->*(pvcInqDimlen[dim]))(dimension_name);
}

// inquire current timestep
int CDI::inq_current_timestep() const {
	return streamInqCurTimestepID(m_file_id) + 1;
}

// inquire time dimension name
void CDI::inq_unlimdim_impl(std::string &result) const {
	result = "time"; // limitation of CDI: cannot set time dimension name
}

// define variable
void CDI::def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
    // No need to define the dimensions as variables
    if (std::find(m_dims_name.begin(), m_dims_name.end(), name) != m_dims_name.end())
    {
        return;
    }
    
    // Define variables
    if (m_vlistID == -1) m_vlistID = vlistCreate();

    if (dims.empty() || dims.size()<2) { // scalar variable
		def_var_scalar_impl(name, nctype, dims);
    } else {                             // multi-dimensional variable
		def_var_multi_impl(name, nctype, dims);
    }
}

void CDI::def_var_scalar_impl(const std::string &name, IO_Type nctype,  const std::vector<std::string> &dims) const {
    if (m_gridsID == -1) {
        m_gridsID = gridCreate(GRID_LONLAT, 1);
        gridDefXsize(m_gridsID, 1);
        gridDefXname(m_gridsID, "x_dummy");
        gridDefYsize(m_gridsID, 1);
		gridDefYname(m_gridsID, "y_dummy");
    }
    if (m_zsID == -1) m_zsID = zaxisCreate(ZAXIS_SURFACE, 1);

    int tsteptype;
    if (dims.empty()) {
        tsteptype = TIME_CONSTANT;
    } else {
	tsteptype = TIME_VARIABLE;
    }

    int varID = vlistDefVar(m_vlistID, m_gridsID, m_zsID, tsteptype);
    vlistDefVarName(m_vlistID, varID, name.c_str());
    int type = pism_type_to_cdi_type(nctype);
    vlistDefVarDatatype(m_vlistID, varID, type);
    m_varsID[name] = varID;
}

void CDI::def_var_multi_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
    int zaxisID = -1;
    int tsteptype = TIME_CONSTANT;

    for (auto d : dims) {
        if        (d == "z") {
            zaxisID = m_zID["z"];
        } else if (d == "zb") {
            zaxisID = m_zbID;
        } else if (d == "time" ) {
            tsteptype = TIME_VARIABLE;
    	}
    }
    if (zaxisID == -1) zaxisID = m_zsID;

    int varID = vlistDefVar(m_vlistID, m_gridID, zaxisID, tsteptype);
    vlistDefVarName(m_vlistID, varID, name.c_str());
    int type = pism_type_to_cdi_type(nctype);
    vlistDefVarDatatype(m_vlistID, varID, type);
    m_varsID[name] = varID;
}

// write spatial dimensions
void CDI::put_dim_x(const std::string &variable_name, const double *op) const {
	gridDefXvals(m_gridID, op);
}

void CDI::put_dim_y(const std::string &variable_name, const double *op) const {
	gridDefYvals(m_gridID, op);
}

void CDI::put_dim_z(const std::string &variable_name, const double *op) const {
	zaxisDefLevels(m_zID[variable_name], op);
}

void CDI::wrapup_put_dim() const {
	pPutDim pFn = &CDI::put_dim_x;
	pvcPutDim.push_back(pFn);
	pFn = &CDI::put_dim_y;
	pvcPutDim.push_back(pFn);
	pFn = &CDI::put_dim_z;
	pvcPutDim.push_back(pFn);
}

// write spatial dimensions and scalars
void CDI::put_vara_double_impl(const std::string &variable_name,
                                  const std::vector<unsigned int> &start,
                                  const std::vector<unsigned int> &count,
                                  const double *op) const {
	// write dimensions values if not done yet
	int dim = m_dimsAxis[variable_name];
	if (!m_gridexist) (this->*(pvcPutDim[dim]))(variable_name, op);
	// if variable_name is a dimension, return
	if (m_dimsAxis.count(variable_name))
		return;
    // write scalar
    int ndims = static_cast<int>(start.size());
    int idxlen=1, j;
    Xt_int *idx;
    for (int j = 0; j < ndims; ++j)
    	idxlen *= count[j];
    idx = (Xt_int*) malloc(idxlen * sizeof(Xt_int));
    for (j=0; j<idxlen; j++)
    	idx[j] = j;
    Xt_idxlist decomp = xt_idxvec_new(idx, idxlen);
    int varid = m_varsID[variable_name];
    size_t nmiss = 0;
    streamWriteVarPart(m_file_id, varid, op, nmiss, decomp);
}

// inquire number of variables
void CDI::inq_nvars_impl(int &result) const {
	result = m_varsID.size();
}

// inquire variable dimensions
void CDI::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
	int varID = m_varsID[variable_name];
	int type = 0;
	int current_grid = vlistInqVarGrid(m_vlistID, varID);
	if (current_grid == m_gridID) {
		int current_z = vlistInqVarZaxis(m_vlistID, varID);
		if (current_z == m_zID["z"]) {
			type = type + 2;
		} else if (current_z == m_zbID) {
			type = type + 4;
		}
                int current_t = vlistInqVarTsteptype(m_vlistID, varID);
		if (current_t == TIME_VARIABLE) {
			type++;
		}
	} else {
		result.clear();
		return;
	}
        // todo: inquire the name - don't write it directly
	switch (type) {
		case 0:
			result.resize(2); result[0] = "y"; result[1] = "x";
			break;
		case 1:
			result.resize(3); result[0] = "time"; result[1] = "y"; result[2] = "x";
			break;
		case 2:
			result.resize(3); result[0] = "z"; result[1] = "y"; result[2] = "x";
			break;
		case 3:
			result.resize(4); result[0] = "time"; result[1] = "z"; result[2] = "y"; result[3] = "x";
			break;
		case 4:
			result.resize(3); result[0] = "zb"; result[1] = "y"; result[2] = "x";
			break;
		case 5:
			result.resize(4); result[0] = "time"; result[1] = "zb"; result[2] = "y"; result[3] = "x";
			break;
	}
}

// inquire variable number of attributes
void CDI::inq_varnatts_impl(const std::string &variable_name, int &result) const {
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	cdiInqNatts(m_vlistID, varID, &result);
}

// inquire variable ID
void CDI::inq_varid_impl(const std::string &variable_name, bool &exists) const {
	exists = m_varsID.count(variable_name) > 0;
}

// inquire variable name
void CDI::inq_varname_impl(unsigned int j, std::string &result) const {
      for (auto &i : m_varsID) {
          if (i.second == j) {
              result = i.first;
              break; // to stop searching
          }
      }
}

// delete variable attribute
void CDI::del_att_impl(const std::string &variable_name, const std::string &att_name) const {
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	cdiDelAtt(m_vlistID, varID, att_name.c_str());
}

// write variable attribute (double)
void CDI::put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type nctype, const std::vector<double> &data) const {
	if (!m_firststep) return;
	// if variable_name is a dimension, return
	if (m_dimsAxis.count(variable_name))
		return;
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	int type = pism_type_to_cdi_type(nctype);
	cdiDefAttFlt(m_vlistID, varID, att_name.c_str(), type, data.size(), &data[0]);
}

// write dimension attribute (text)
void CDI::put_att_text_units_x_impl(const std::string &variable_name, const std::string &value) const {
	gridDefXunits(m_gridID, value.c_str());
}

void CDI::put_att_text_longname_x_impl(const std::string &variable_name, const std::string &value) const {
	gridDefXlongname(m_gridID, value.c_str());
}

void CDI::put_att_text_units_y_impl(const std::string &variable_name, const std::string &value) const {
	gridDefYunits(m_gridID, value.c_str());
}

void CDI::put_att_text_longname_y_impl(const std::string &variable_name, const std::string &value) const {
	gridDefYlongname(m_gridID, value.c_str());
}

void CDI::put_att_text_units_z_impl(const std::string &variable_name, const std::string &value) const {
	zaxisDefUnits(m_zID[variable_name], value.c_str());
}

void CDI::put_att_text_longname_z_impl(const std::string &variable_name, const std::string &value) const {
	zaxisDefLongname(m_zID[variable_name], value.c_str());
}

void CDI::wrapup_put_att_text() const {
	m_DimAtt["units"] = 0;
	m_DimAtt["long_name"] = 1;
	pvcPutAttT[0][0] = &CDI::put_att_text_units_x_impl;
	pvcPutAttT[0][1] = &CDI::put_att_text_longname_x_impl;
	pvcPutAttT[1][0] = &CDI::put_att_text_units_y_impl;
	pvcPutAttT[1][1] = &CDI::put_att_text_longname_y_impl;
	pvcPutAttT[2][0] = &CDI::put_att_text_units_z_impl;
	pvcPutAttT[2][1] = &CDI::put_att_text_longname_z_impl;
}

void CDI::put_att_text_impl(const std::string &variable_name,
                                const std::string &att_name,
                                const std::string &value) const {
	if (!m_firststep) return;
	// basic checks
	if (value.empty() || att_name.empty()) return;
	// write dimension attribute
    if (m_dimsAxis.count(variable_name)) {
    	int dim = m_dimsAxis[variable_name];
    	int att = m_DimAtt[att_name];
    	(this->*(pvcPutAttT[dim][att]))(variable_name, value);
        return;
    }

    // write variable attribute
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	cdiDefAttTxt(m_vlistID, varID, att_name.c_str(), value.size(), value.c_str());
}

// inquire attribute type
void CDI::inq_atttype_impl(const std::string &variable_name,
                               const std::string &att_name,
                               IO_Type &result) const {
	// find variable ID
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}

	// find attribute type
	int natt, n, atype, alen, cditype;
	char name[CDI_MAX_NAME];
	std::string aname;
	if (att_name=="history" && varID==CDI_GLOBAL) {
		cditype = CDI_DATATYPE_TXT;
	} else {
		cdiInqNatts(m_vlistID, varID, &natt);
		for (n=0; n<natt; n++) {
			inq_att_impl(varID, n, name, &atype, &alen);
			if (aname==att_name) {
				cditype = atype;
			}
		}
	}
	result = cdi_type_to_pism_type(cditype);
}

// inquire attribute name
void CDI::inq_attname_impl(const std::string &variable_name,
                               unsigned int n,
                               std::string &result) const {
	char name[CDI_MAX_NAME];
	int atype, alen;
	// find variable ID
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	// find attribute name
	inq_att_impl(varID, n, name, &atype, &alen);
	result = name;
}

void CDI::inq_att_impl(int varID, int attnum, char* attname, int *atttype, int *attlen) const {
	cdiInqAtt(m_vlistID, varID, attnum, attname, atttype, attlen);
}

// get variable attribute (double)
void CDI::get_att_double_impl(const std::string &variable_name,
                                  const std::string &att_name,
                                  std::vector<double> &result) const {
	// find variable ID
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}

	// find attribute length
	int natt, n, atype, alen, cdilen;
	char name[CDI_MAX_NAME];
	std::string aname;
	cdiInqNatts(m_vlistID, varID, &natt);
	for (n=0; n<natt; n++) {
		inq_att_impl(varID, n, name, &atype, &alen);
		aname = name;
		if (aname==att_name) {
			cdilen = alen;
		}
	}

	if (cdilen == 0) {
    	result.clear();
    	return;
  	}
  	result.resize(cdilen);
	
	// read attribute
	cdiInqAttFlt(m_vlistID, varID, att_name.c_str(), cdilen, std::addressof(result[0]));
}

// get variable attribute (text)
void CDI::get_att_text_impl(const std::string &variable_name, 
						    const std::string &att_name, 
						    std::string &result) const {
	// find variable ID
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}

	// find attribute length
	int natt, n, atype, alen, cdilen;
	char name[CDI_MAX_NAME];
	std::string aname;
	if (att_name=="history" && varID==CDI_GLOBAL) {
		cdilen = streamInqHistorySize(m_file_id);
	} else {
		cdiInqNatts(m_vlistID, varID, &natt);
		for (n=0; n<natt; n++) {
			inq_att_impl(varID, n, name, &atype, &alen);
			if (aname==att_name) {
				cdilen = alen;
			}
		}
	}

	if (cdilen == 0) {
    	result.clear();
    	return;
  	}
  	result.resize(cdilen);

  	// read attribute
  	if (att_name=="history" && varID==CDI_GLOBAL) {
  		streamInqHistoryString(m_file_id, std::addressof(result[0]));
  	} else {
  		cdiInqAttTxt(m_vlistID, varID, att_name.c_str(), cdilen, std::addressof(result[0]));
  	}
}

// create main grid
void CDI::create_grid_impl(int lengthx, int lengthy) const {
	if (m_gridID ==-1) m_gridID = gridCreate(GRID_LONLAT, lengthx*lengthy);
}

// define timestep
void CDI::define_timestep_impl(int tsID) const {
	streamDefTimestep(m_file_id, tsID);
	if (tsID==0) m_firststep = false;
}

// write variables
void CDI::write_darray_impl(const std::string &variable_name,
                                   const IceGrid &grid,
                                   unsigned int z_count,
                                   unsigned int record,
                                   const double *input) {
	// transpose input data
	int dim = grid.local_length((int)z_count);
	double *inputIO;
	inputIO = (double*) malloc(dim * sizeof(double));
	grid.io_transpose(input, inputIO, (int)z_count);

	int varid = m_varsID[variable_name];
	
	// create decomposition if new
	Xt_idxlist decompid = grid.yaxt_decomposition((int)z_count);
	size_t nmiss = 0;
	if (!m_beforediag || m_diagvars.count(variable_name)==0) {
		streamWriteVarPart(m_file_id, varid, inputIO, nmiss, decompid);
	}
}

std::map<std::string, int> CDI::get_var_map_impl() {
	return m_varsID;
}

void CDI::def_vlist_impl() const {
	if (m_firststep) streamDefVlist(m_file_id, m_vlistID);
}

void CDI::set_diagvars_impl(const std::set<std::string> &variables) const {
	m_diagvars = variables;
}

void CDI::set_bdiag_impl(bool value) const {
	m_beforediag = value;
}

int CDI::get_ncstreamID_impl() const {
	return m_file_id;
}

int CDI::get_ncvlistID_impl() const {
	return m_vlistID;
}

// Not used
void CDI::sync_impl() const {
}

void CDI::enddef_impl() const {
}

void CDI::redef_impl() const {
}

void CDI::get_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 double *ip) const {
        (void) variable_name;
        (void) start;
        (void) count;
        (void) ip;
}

void CDI::get_varm_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap,
                                 double *ip) const {
        (void) variable_name;
        (void) start;
        (void) count;
        (void) imap;
        (void) ip;
}



void CDI::set_fill_impl(int fillmode, int &old_modep) const {
        (void) fillmode;
        (void) old_modep;
}

} // end of namespace io
} // end of namespace pism
