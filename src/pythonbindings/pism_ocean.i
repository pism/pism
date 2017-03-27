%{
#include "coupler/ocean/POConstant.hh"
%}

%shared_ptr(pism::Diag<pism::ocean::OceanModel>)
%template(OceanDiagnostic) pism::Diag<pism::ocean::OceanModel>;
%shared_ptr(pism::ocean::OceanModel)
%shared_ptr(pism::ocean::PO_sea_level)
%shared_ptr(pism::ocean::PO_shelf_base_temperature)
%shared_ptr(pism::ocean::PO_shelf_base_mass_flux)
%shared_ptr(pism::ocean::PO_melange_back_pressure_fraction)
%include "coupler/PISMOcean.hh"

%shared_ptr(pism::ocean::Constant)
%rename(OceanConstant) pism::ocean::Constant;
%include "coupler/ocean/POConstant.hh"
