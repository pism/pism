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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <mpi.h>
#include <string.h>

#include "CDI.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/io/pism_cdi_type_conversion.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/error_handling.hh"
#include "pism/external/calcalcs/calcalcs.h"
#include "pism/util/cdipio/Idxlist.hh"

extern "C" {
#include "cdi.h"
#include "cdipio.h"
#include "yaxt.h"
}

namespace pism {
namespace io {

CDI::CDI(MPI_Comm c) : NCFile(c) {
  m_beforediag = true;
}

CDI::~CDI() {
  // empty
}

void CDI::open_impl(const std::string &fname,
                    IO_Mode mode,
                    int FileID,
                    const std::map<std::string, AxisType> &dimsa) {
  // FIXME: in general the assumption below is incorrect
  //
  // the file is already created and opened - restore file info into the class
  if (mode == PISM_READONLY) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "file reading not supported with CDI-PIO in PISM");
  }
  m_file_id = FileID;
  m_vlistID = streamInqVlist(m_file_id);
  m_tID     = vlistInqTaxis(m_vlistID);
  map_varsID();
  map_zaxisID();
  m_gridID   = vlistGrid(m_vlistID, 0);
  m_dimsAxis = dimsa;
}

void CDI::map_varsID() const {
  int Nvars, varID;
  char VarName[CDI_MAX_NAME];
  Nvars = vlistNvars(m_vlistID);
  for (varID = 0; varID < Nvars; varID++) {
    vlistInqVarName(m_vlistID, varID, VarName);
    m_varsID[VarName] = varID;
  }
}

void CDI::map_zaxisID() const {
  int zaxisID;
  char zname[CDI_MAX_NAME];
  // find number of zaxis
  int nz = vlistNzaxis(m_vlistID);
  // find zaxisID and zaxisName
  for (int n = 0; n < nz; n++) {
    zaxisID = vlistZaxis(m_vlistID, n);
    zaxisInqName(zaxisID, zname);
    m_zID[zname] = zaxisID;
  }
  m_zsID = m_zID["zs"];
}

void CDI::create_impl(const std::string &filename, int FileID, const std::string &filetype) {
  // FIXME: parameter FileID is not used
  int itype = define_filetype(filetype);
  m_file_id = streamOpenWrite(filename.c_str(), itype);
  m_tID     = -1;
  m_zsID    = -1;
  m_gridID  = -1;
  m_vlistID = -1;
}

int CDI::define_filetype(const std::string &filetype) {
  std::map<std::string, int> types =
    {
     {"CDI_FILETYPE_NC", CDI_FILETYPE_NC},
     {"CDI_FILETYPE_NC2", CDI_FILETYPE_NC2},
     {"CDI_FILETYPE_NC4", CDI_FILETYPE_NC4},
     {"CDI_FILETYPE_NC4C", CDI_FILETYPE_NC4C}
    };

  if (types.find(filetype) == types.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid CDI-PIO file type %s", filetype.c_str());
  }

  return types[filetype];
}

void CDI::close_impl() {
  m_file_id = -1;
  m_varsID.clear();
  m_dimsAxis.clear();
  m_zID.clear();
  m_diagvars.clear();
}

void CDI::def_vlist() const {
  if (m_vlistID == -1) {
    // create variable list
    m_vlistID = vlistCreate();
    // create dummy grid to handle scalar values
    m_gridsID = gridCreate(GRID_LONLAT, 1);
    gridDefXsize(m_gridsID, 1);
    gridDefXname(m_gridsID, "x_dummy");
    gridDefYsize(m_gridsID, 1);
    gridDefYname(m_gridsID, "y_dummy");
  }
}

void CDI::def_zs() const {
  // create surface Z axis (only if it was not done before)
  if (m_zsID == -1) {
    m_zsID      = zaxisCreate(ZAXIS_SURFACE, 1);
    m_zID["zs"] = m_zsID;
    zaxisDefName(m_zsID, "zs");
  }
}

void CDI::def_dim_impl(const std::string &name, size_t length, AxisType dim) const {

  def_vlist();
  def_zs();

  switch (dim) {
  case X_AXIS:
    {
      gridDefXsize(m_gridID, (int)length);
      gridDefXname(m_gridID, name.c_str());
      break;
    }
  case Y_AXIS:
    {
      gridDefYsize(m_gridID, (int)length);
      gridDefYname(m_gridID, name.c_str());
      break;
    }
  case Z_AXIS:
    {
      // define z axis only if it's new
      if (!m_zID.count(name)) {
        m_zID[name] = zaxisCreate(ZAXIS_GENERIC, (int)length);
        zaxisDefName(m_zID[name], name.c_str());
      }
      break;
    }
  case T_AXIS:
    {
      // define time axis if it was not done before
      if (m_tID == -1) {
        m_tID = taxisCreate(TAXIS_ABSOLUTE);
        taxisDefCalendar(m_tID, m_cdi_calendar);
        vlistDefTaxis(m_vlistID, m_tID);
      }
      break;
    }
  default:
    {
      // define grid for 1D data
      m_gridgID = gridCreate(GRID_LONLAT, length);
      gridDefXsize(m_gridgID, length);
      gridDefXname(m_gridgID, name.c_str());
      gridDefYsize(m_gridgID, 1);
      gridDefYname(m_gridgID, "y_dummy");
      break;
    }
  }
  m_dimsAxis[name] = dim;
}

// set file calendar (CDIPIO needs to know calendar type)
void CDI::set_calendar_impl(double year_length, const std::string &calendar_string) const {
  m_year_length = year_length;
  if (calendar_string == "gregorian" || calendar_string == "standard") {
    m_days_year    = m_year_length / 86400;
    m_cdi_calendar = CALENDAR_STANDARD;
  } else if (calendar_string == "proleptic_gregorian") {
    m_days_year    = m_year_length / 86400;
    m_cdi_calendar = CALENDAR_PROLEPTIC;
  } else if (calendar_string == "365_day" || calendar_string == "noleap") {
    m_days_year    = 365;
    m_cdi_calendar = CALENDAR_365DAYS;
  } else if (calendar_string == "360_day") {
    m_days_year    = 360;
    m_cdi_calendar = CALENDAR_360DAYS;
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "unsupported calendar: %s",
                                  calendar_string.c_str());
  }
  m_calendar_string = calendar_string;
}

// evaluate year from time input
double CDI::year_calendar(double time) const {
  return abs(time / m_year_length);
}

// evaluate month and day
void CDI::monthday_calendar(int year, int doy, int *month, int *day) const {
  calcalcs_cal *calendar = ccs_init_calendar(m_calendar_string.c_str());
  assert(calendar != NULL);
  int cal = ccs_doy2date(calendar, year, doy, month, day);
  assert(cal == 0);
  ccs_free_calendar(calendar);
}

// evaluate hours, minutes and seconds
long int CDI::day_calendar(double nyearsf) const {
  long int seconds = round(nyearsf * 86400);
  long int minutes, hours;
  hours   = seconds / 3600;
  minutes = (seconds - (3600 * hours)) / 60;
  seconds = (seconds - (3600 * hours) - minutes * 60);
  return hours * 10000 + minutes * 100 + seconds;
}

// define time (CDI support specific time format)
void CDI::def_ref_date_impl(double time) const {
  // conversion of "time" input into CDI format 
  int month = 0, day = 0;
  double nyearsf = year_calendar(time);
  int doy        = int((nyearsf - (long int)nyearsf) * m_days_year);
  double dayf    = ((nyearsf - (long int)nyearsf) * m_days_year) - doy;
  if (doy != 0) {
    monthday_calendar(int(nyearsf), doy, &month, &day);
  }
  // define time in CDI format
  int ref_date = (int)nyearsf * 10000 + month * 100 + day;
  int sgn      = time >= 0 ? 1 : -1;
  taxisDefVdate(m_tID, sgn * ref_date);
  long int daytime = day_calendar(dayf);
  taxisDefVtime(m_tID, daytime);
}

// inquire if a dimension exists
void CDI::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
  exists = m_dimsAxis.count(dimension_name) > 0;
}

// inquire dimension length
void CDI::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
  int dim = m_dimsAxis[dimension_name];
  switch (dim) {
  case X_AXIS:
    result = gridInqXsize(m_gridID);
    break;
  case Y_AXIS:
    result = gridInqYsize(m_gridID);
    break;
  case Z_AXIS:
    result = zaxisInqSize(m_zID[dimension_name]);
    break;
  case T_AXIS:
    result = streamInqCurTimestepID(m_file_id) + 1 + 1;
    break;
  default:
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid dim == %d", dim);
  }
}

// inquire time dimension name
void CDI::inq_unlimdim_impl(std::string &result) const {
  result = "time"; // limitation of CDI: cannot set time dimension name
}

// define variable
void CDI::def_var_impl(const std::string &name, IO_Type nctype,
                       const std::vector<std::string> &dims) const {
  // No need to define dimensions as variables (CDI does not support it)
  if (m_dimsAxis.count(name)) {
    return;
  }

  // Define variables
  def_vlist();
  int tdim = 0;
  if (std::find(dims.begin(), dims.end(), "time") != dims.end()) {
    tdim = 1;
  }

  if (dims.empty() || dims.size() - tdim == 0) { // scalar variable
    def_var_scalar(name, nctype, dims);
  } else if (dims.size() - tdim == 1) { // time-dependent scalar variable
    def_var_mscalar(name, nctype, dims);
  } else { // multi-dimensional variable
    def_var_multi(name, nctype, dims);
  }
}

void CDI::def_var_scalar(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
  // define surface Z axis (if not done before)
  def_zs();

  // inquire if variable is time dependent or not
  int tsteptype = TIME_CONSTANT;
  for (auto d : dims) {
    auto it = m_dimsAxis.find(d);
    if (it != m_dimsAxis.end() && it->second == T_AXIS) {
      tsteptype = TIME_VARIABLE;
    }
  }

  // define variable - dummy scalar grid is used
  int varID = vlistDefVar(m_vlistID, m_gridsID, m_zsID, tsteptype);
  vlistDefVarName(m_vlistID, varID, name.c_str());
  int type = pism_type_to_cdi_type(nctype);
  vlistDefVarDatatype(m_vlistID, varID, type);
  // add variable ID to the map - variable ID can not be inquired (needs to be saved)
  m_varsID[name] = varID;
}

void CDI::def_var_mscalar(const std::string &name, IO_Type nctype,
                          const std::vector<std::string> &dims) const {
  // define surface Z axis (if not done before)
  def_zs();

  // inquire if variable is time dependent or not
  int tsteptype = TIME_CONSTANT;
  for (auto d : dims) {
    auto it = m_dimsAxis.find(d);
    if (it != m_dimsAxis.end() && it->second == T_AXIS) {
      tsteptype = TIME_VARIABLE;
    }
  }

  // define variable - dummy 1D grid is used
  int varID = vlistDefVar(m_vlistID, m_gridgID, m_zsID, tsteptype);
  vlistDefVarName(m_vlistID, varID, name.c_str());
  int type = pism_type_to_cdi_type(nctype);
  vlistDefVarDatatype(m_vlistID, varID, type);
  // add variable ID to the map - variable ID can not be inquired (needs to be saved)
  m_varsID[name] = varID;
}


void CDI::def_var_multi(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
  int zaxisID   = -1;
  int tsteptype = TIME_CONSTANT;

  // inquire if variable is time dependent or not and the associated Z axis
  for (auto d : dims) {
    if (m_dimsAxis[d] == T_AXIS) {
      tsteptype = TIME_VARIABLE;
    }
    if (m_dimsAxis[d] == Z_AXIS) {
      if (m_zID.count(d)) {
        zaxisID = m_zID[d];
      }
    }
  }

  if (zaxisID == -1) {
    zaxisID = m_zsID;
  }

  // define variable
  int varID = vlistDefVar(m_vlistID, m_gridID, zaxisID, tsteptype);
  vlistDefVarName(m_vlistID, varID, name.c_str());
  int type = pism_type_to_cdi_type(nctype);
  vlistDefVarDatatype(m_vlistID, varID, type);
  // add variable ID to the map - variable ID can not be inquired (needs to be saved)
  m_varsID[name] = varID;
}

// write spatial dimensions and scalars
void CDI::put_vara_double_impl(const std::string &variable_name, const std::vector<unsigned int> &start,
                               const std::vector<unsigned int> &count, const double *op) const {
  // write dimensions values if not done yet
  if (m_dimsAxis.count(variable_name)) {
    int dim = m_dimsAxis[variable_name];
    switch (dim) {
    case X_AXIS:
      gridDefXvals(m_gridID, op);
      break;
    case Y_AXIS:
      gridDefYvals(m_gridID, op);
      break;
    case Z_AXIS:
    default:
      zaxisDefLevels(m_zID[variable_name], op);
      break;
    }
    return;
  }

  // write scalar
  int ndims  = static_cast<int>(start.size());
  // define dummy YAXT decomp to write scalar variables with CDIPIO
  int idxlen = 1;
  Xt_int *idx;
  for (int j = 0; j < ndims; ++j)
    idxlen *= count[j];
  idx = (Xt_int *)malloc(idxlen * sizeof(Xt_int));
  for (int j = 0; j < idxlen; j++)
    idx[j] = j;
  Xt_idxlist decomp = xt_idxvec_new(idx, idxlen);
  free(idx);
  // inquire variable ID from the map
  int varid         = m_varsID[variable_name];
  size_t nmiss      = 0;
  // write scalar variable
  streamWriteVarPart(m_file_id, varid, op, nmiss, decomp);
}

// inquire number of variables
void CDI::inq_nvars_impl(int &result) const {
  result = m_varsID.size();
}

// inquire variable dimensions
void CDI::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
  std::string name;
  int varID        = m_varsID[variable_name];
  int current_grid = vlistInqVarGrid(m_vlistID, varID);
  char xname[CDI_MAX_NAME];

  if (current_grid == m_gridID) {
    // insert x dim
    auto it = result.begin();
    gridInqXname(m_gridID, xname);
    name = xname;
    result.insert(it, name);
    // insert y dim
    it = result.begin();
    gridInqYname(m_gridID, xname);
    name = xname;
    result.insert(it, name);
    // insert z dim
    name.clear();
    int current_z = vlistInqVarZaxis(m_vlistID, varID);
    for (auto &i : m_zID) {
      if (i.second == current_z) {
        name = i.first;
        break; // to stop searching
      }
    }
    if (!name.empty()) {
      it = result.begin();
      result.insert(it, name);
    }
    // insert time dim
    if (vlistInqVarTsteptype(m_vlistID, varID) == TIME_VARIABLE) {
      it = result.begin();
      result.insert(it, "time");
    }
  } else {
    result.clear();
    return;
  }
}

// inquire variable ID using map
int CDI::var_id(const std::string &name) const {
  if (name == "PISM_GLOBAL") {
    return CDI_GLOBAL;
  } else {
    return m_varsID[name];
  }
}

// inquire variable number of attributes
void CDI::inq_varnatts_impl(const std::string &variable_name, int &result) const {
  int varID = var_id(variable_name);
  cdiInqNatts(m_vlistID, varID, &result);
}

// inquire variable ID
void CDI::inq_varid_impl(const std::string &variable_name, bool &exists) const {
  exists = m_varsID.count(variable_name) > 0;
}

// inquire variable name
void CDI::inq_varname_impl(unsigned int j, std::string &result) const {
  for (auto &i : m_varsID) {
    if (i.second == (int)j) {
      result = i.first;
      break; // to stop searching
    }
  }
}

// delete variable attribute
void CDI::del_att_impl(const std::string &variable_name, const std::string &att_name) const {
  int varID = var_id(variable_name);
  cdiDelAtt(m_vlistID, varID, att_name.c_str());
}

// write variable attribute (double)
void CDI::put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type nctype,
                              const std::vector<double> &data) const {
  // if variable_name is a dimension, return
  if (m_dimsAxis.count(variable_name)) {
    return;
  }

  int varID = var_id(variable_name);
  int type = pism_type_to_cdi_type(nctype);
  cdiDefAttFlt(m_vlistID, varID, att_name.c_str(), type, data.size(), &data[0]);
}

// write variable attribute (text)
void CDI::put_att_text_impl(const std::string &variable_name, const std::string &att_name,
                            const std::string &value) const {
  // skip empty attributes
  if (value.empty() or att_name.empty()) {
    return;
  }

  bool dimension = m_dimsAxis.count(variable_name) > 0;
  // write dimension attribute (limited supported dimension attributes)
  if (dimension) {
    int type = m_dimsAxis[variable_name];
    switch (type) {
    case X_AXIS:
      if (att_name == "units") {
        gridDefXunits(m_gridID, value.c_str());
      }
      if (att_name == "long_name") {
        gridDefXlongname(m_gridID, value.c_str());
      }
      break;
    case Y_AXIS:
      if (att_name == "units") {
        gridDefYunits(m_gridID, value.c_str());
      }
      if (att_name == "long_name") {
        gridDefYlongname(m_gridID, value.c_str());
      }
      break;
    case Z_AXIS:
      if (att_name == "units") {
        zaxisDefUnits(m_zID[variable_name], value.c_str());
      }
      if (att_name == "long_name") {
        zaxisDefLongname(m_zID[variable_name], value.c_str());
      }
      break;
    default:
      break;
    }
    // all other dimension attributes are ignored
    return;
  }


  // write variable attribute
  int varid = var_id(variable_name);

  cdiDefAttTxt(m_vlistID, varid, att_name.c_str(), value.size(), value.c_str());
}

// inquire attribute type
void CDI::inq_atttype_impl(const std::string &variable_name, const std::string &att_name, IO_Type &result) const {
  // find variable ID
  int varID = var_id(variable_name);

  // find attribute type
  int natt, n, atype, alen, cditype;
  char name[CDI_MAX_NAME];
  std::string aname;
  if (att_name == "history" && varID == CDI_GLOBAL) {
    cditype = CDI_DATATYPE_TXT;
  } else {
    cdiInqNatts(m_vlistID, varID, &natt);
    for (n = 0; n < natt; n++) {
      inq_att_impl(varID, n, name, &atype, &alen);
      if (aname == att_name) {
        cditype = atype;
      }
    }
  }
  result = cdi_type_to_pism_type(cditype);
}

// inquire attribute name
void CDI::inq_attname_impl(const std::string &variable_name, unsigned int n, std::string &result) const {
  char name[CDI_MAX_NAME];
  int atype, alen;
  // find variable ID
  int varID = var_id(variable_name);
  // find attribute name
  inq_att_impl(varID, n, name, &atype, &alen);
  result = name;
}

// inquire attribute name, time and length
void CDI::inq_att_impl(int varID, int attnum, char *attname, int *atttype, int *attlen) const {
  cdiInqAtt(m_vlistID, varID, attnum, attname, atttype, attlen);
}

// get variable attribute (double)
void CDI::get_att_double_impl(const std::string &variable_name, const std::string &att_name,
                              std::vector<double> &result) const {
  // find variable ID
  int varID = var_id(variable_name);

  // find attribute length
  int natt, n, atype, alen, cdilen;
  char name[CDI_MAX_NAME];
  std::string aname;
  cdiInqNatts(m_vlistID, varID, &natt);
  for (n = 0; n < natt; n++) {
    inq_att_impl(varID, n, name, &atype, &alen);
    aname = name;
    if (aname == att_name) {
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
void CDI::get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const {
  // find variable ID
  int varID = var_id(variable_name);

  // find attribute length
  int natt, n, atype, alen, cdilen;
  char name[CDI_MAX_NAME];
  std::string aname;
  if (att_name == "history" && varID == CDI_GLOBAL) {
    cdilen = streamInqHistorySize(m_file_id);
  } else {
    cdiInqNatts(m_vlistID, varID, &natt);
    for (n = 0; n < natt; n++) {
      inq_att_impl(varID, n, name, &atype, &alen);
      if (aname == att_name) {
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
  if (att_name == "history" && varID == CDI_GLOBAL) {
    streamInqHistoryString(m_file_id, std::addressof(result[0]));
  } else {
    cdiInqAttTxt(m_vlistID, varID, att_name.c_str(), cdilen, std::addressof(result[0]));
  }
}

// create main grid
void CDI::create_grid_impl(int lengthx, int lengthy) const {
  if (m_gridID == -1) {
    // We should use GRID_PROJECTION, but CDI halts with the error "unimplemented grid
    // type: 12" in cdiPioQueryVarDims().
    m_gridID = gridCreate(GRID_LONLAT, lengthx * lengthy);
  }
  // FIXME: that happens if this is called twice, but with different lengthx and lengthy?
}

// define timestep
void CDI::define_timestep_impl(int tsID) const {
  streamDefTimestep(m_file_id, tsID);
}

// write variables
void CDI::write_darray_impl(const std::string &variable_name,
                            const IceGrid &grid, unsigned int z_count,
                            unsigned int record, const double *input) {
  (void) record;
  // CDI cannot write an arbitrary record in the file, so "record" is ignored

  // transpose input data (limitation of CDIPIO)
  int dim = grid.local_length(z_count);
  std::vector<double> buffer(dim);
  grid.io_transpose(input, buffer.data(), (int)z_count);

  int varid = var_id(variable_name);
  assert(varid != CDI_GLOBAL);

  // create decomposition if new
  Xt_idxlist decompid = *grid.yaxt_decomposition(z_count);
  size_t nmiss        = 0;
  // write variable if it is written once or if this is the last call
  // CDIPIO does not support writing the same variable multiple times
  // at the same time step
  if (!m_beforediag || m_diagvars.count(variable_name) == 0) {
    streamWriteVarPart(m_file_id, varid, buffer.data(), nmiss, decompid);
  }
}

std::map<std::string, int> CDI::get_var_map_impl() {
  return m_varsID;
}

std::map<std::string, AxisType> CDI::get_dim_map_impl() {
  return m_dimsAxis;
}

// define variables list (if not done before)
void CDI::def_vlist_impl() const {
  if (streamInqVlist(m_file_id) == -1) {
    streamDefVlist(m_file_id, m_vlistID);
  }
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

void CDI::get_vara_double_impl(const std::string &variable_name, const std::vector<unsigned int> &start,
                               const std::vector<unsigned int> &count, double *ip) const {
  (void)variable_name;
  (void)start;
  (void)count;
  (void)ip;
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "file reading not supported with CDI-PIO in PISM");
}

void CDI::get_varm_double_impl(const std::string &variable_name, const std::vector<unsigned int> &start,
                               const std::vector<unsigned int> &count, const std::vector<unsigned int> &imap,
                               double *ip) const {
  (void)variable_name;
  (void)start;
  (void)count;
  (void)imap;
  (void)ip;
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "file reading not supported with CDI-PIO in PISM");
}

void CDI::set_fill_impl(int fillmode, int &old_modep) const {
  (void)fillmode;
  (void)old_modep;
}

} // end of namespace io
} // end of namespace pism
