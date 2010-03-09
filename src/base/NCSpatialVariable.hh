
#ifndef __NCSpatialVariable
#define __NCSpatialVariable

#include "NCVariable.hh"
#include "LocalInterpCtx.hh"

//! Spatial NetCDF variable (corresponding to a 2D or 3D scalar field).
class NCSpatialVariable : public NCVariable {
public:
  NCSpatialVariable();
  virtual void init(string name, IceGrid &g, GridType d);
  virtual PetscErrorCode read(const char filename[], unsigned int time, Vec v);
  virtual PetscErrorCode reset();
  virtual PetscErrorCode write(const char filename[], nc_type nctype,
			       bool write_in_glaciological_units, Vec v);
  virtual PetscErrorCode regrid(const char filename[], LocalInterpCtx &lic,
				bool critical, bool set_default_value,
				PetscScalar default_value, MaskInterp *, Vec v);
  virtual PetscErrorCode to_glaciological_units(Vec v);
  bool time_independent;
protected:
  GridType dims;
  IceGrid *grid;
  PetscErrorCode define(int ncid, nc_type nctype, int &varid);
  PetscErrorCode report_range(Vec v, bool found_by_standard_name);
  PetscErrorCode change_units(Vec v, utUnit *from, utUnit *to);
  PetscErrorCode check_range(Vec v);
};

#endif	// __NCSpatialVariable
