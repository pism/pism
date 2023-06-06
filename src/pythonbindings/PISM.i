// Copyright (C) 2011--2023 David Maxwell and Constantine Khroulev
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

%module(directors="1") cpp
%feature("autodoc", "2");

/* Don't warn about
 * - 312,325 Nested classes not currently supported (ignored).
 * - 503. Can't wrap 'identifier' unless renamed to a valid identifier.
 * - 512. Overloaded declaration const ignored. Non-const method at file:line used.
 */
#pragma SWIG nowarn=312,325,503,512

%{
// The material in this section is included verbatim in the C++ source code generated by SWIG.
// The necessary header files required to compile must be included.
// This list is NOT the whole set of headers being wrapped; it is just the list of includes that
// draws in all the other needed includes as well. See the end of this file for the list
// of PISM headers being wrapped.

#include "util/interpolation.hh"

#include "util/pism_utilities.hh"

#include "util/Units.hh"
#include "pism_python.hh"

#include "geometry/grounded_cell_fraction.hh"
#include "util/cell_type.hh"
#include "basalstrength/basal_resistance.hh"
#include "basalstrength/MohrCoulombYieldStress.hh"
#include "util/error_handling.hh"
#include "util/Diagnostic.hh"
#include "util/Config.hh"

#if (Pism_USE_JANSSON==1)
#include "util/ConfigJSON.hh"
#endif

#include "util/MaxTimestep.hh"
#include "stressbalance/timestepping.hh"
#include "util/Context.hh"
#include "util/Profiling.hh"

#include "util/projection.hh"
#include "energy/bootstrapping.hh"
#include "util/node_types.hh"

#include "util/label_components.hh"
%}

// Tell SWIG that the following variables are truly constant
%immutable pism::revision;
%immutable pism::config_file;
%immutable pism::petsc_configure_flags;
%immutable pism::petsc4py_version;
%immutable pism::swig_version;
%immutable pism::cmake_version;
%include "pism/pism_config.hh"

// Include petsc4py.i so that we get support for automatic handling of PetscErrorCode return values
%include "petsc4py/petsc4py.i"

%include "pism_exception.i"

// Automatic conversions between std::string and python string arguments and return values
%include std_string.i
// Conversions between python lists and certain STL vectors
%include std_array.i
%include std_vector.i
%include std_set.i
%include std_map.i

%include <std_shared_ptr.i>

// Add a PISM class
%define pism_class(name, header)
%{
  #include header
%}
%shared_ptr(name)
%include header
%enddef

%template(SizetVector) std::vector<size_t>;
%template(IntVector) std::vector<int>;
%template(UnsignedIntVector) std::vector<unsigned int>;
%template(DoubleVector) std::vector<double>;
%template(StringVector) std::vector<std::string>;
%template(StringSet) std::set<std::string>;
%template(DoubleVectorMap) std::map<std::string, std::vector<double> >;
%template(BoolMap) std::map<std::string, bool >;
%template(StringMap) std::map<std::string, std::string>;
%template(DiagnosticMap) std::map<std::string, std::shared_ptr<pism::Diagnostic> >;
%template(SizeDoubleMap) std::map<size_t, double>;

// Why did I include this?
%include "cstring.i"

/* Type map for treating reference arguments as output. */
%define %Pism_reference_output_typemaps(TYPE)
%typemap(in, numinputs=0,noblock=1) TYPE & OUTPUT (TYPE temp) {
    $1 = &temp;
}
%typemap(argout,noblock=1) TYPE & OUTPUT
{
    %append_output(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, 0 | %newpointer_flags));
};
%enddef

/* Tell SWIG that reference arguments are always output. */
%define %Pism_reference_is_always_output(TYPE)
%Pism_reference_output_typemaps(TYPE);
%apply TYPE & OUTPUT { TYPE &}
%enddef

%typemap(in, numinputs=0,noblock=1) bool & OUTPUT (bool temp = false) {
    $1 = &temp;
}

%typemap(argout,noblock=1) bool & OUTPUT
{
    %append_output(SWIG_From(bool)(*$1));
};

%typemap(in, numinputs=0,noblock=1) std::string& result (std::string temp) {
    $1 = &temp;
}

%typemap(in, numinputs=0,noblock=1) std::string& OUTPUT (std::string temp) {
    $1 = &temp;
}

%typemap(argout,noblock=1) std::string & OUTPUT
{
    %append_output(SWIG_FromCharPtr((*$1).c_str()));
}

%apply std::string &OUTPUT { std::string &result}

%typemap(in, numinputs=0,noblock=1) std::vector<double> & OUTPUT (std::vector<double> temp) {
    $1 = &temp;
}

%typemap(argout,noblock=1) std::vector<double> & OUTPUT
{
    int len;
    len = $1->size();
    $result = PyList_New(len);
     int i;
     for (i=0; i<len; i++) {
         PyList_SetItem($result, i, PyFloat_FromDouble((*$1)[i]));
     }
}

%apply std::vector<double> & OUTPUT {std::vector<double> &result};
%apply std::vector<std::string> & OUTPUT {std::vector<std::string> & result};

%apply int &OUTPUT {int &result};
%apply int *OUTPUT {int *out_mask};

%apply double & OUTPUT {double & result};
%apply double & OUTPUT {double & out};
%apply double * OUTPUT {double * result};
%apply bool & OUTPUT {bool & is_set, bool & result, bool & flag, bool & success};

// The SWIG built-in typecheck for a const char [] (used, e.g., with overloaded methods) checks that
// the string is zero length. So we have this bug fix from SWIG developer William Fulton here.
%typemap(typecheck,noblock=1,precedence=SWIG_TYPECHECK_STRING, fragment="SWIG_AsCharPtrAndSize") const char[] {
 int res = SWIG_AsCharPtrAndSize($input, 0, NULL, 0);
 $1 = SWIG_CheckState(res);
}


/* PISM header with no dependence on other PISM headers. */
%ignore print_vector;
%include "util/pism_utilities.hh"
%include "util/interpolation.hh"

%shared_ptr(pism::StringLogger);
pism_class(pism::Logger, "pism/util/Logger.hh");

%include pism_options.i

%ignore pism::Vector2d::operator=;
%include "util/Vector2d.hh"

%ignore pism::units::Unit::operator=;
%rename(UnitSystem) pism::units::System;
%rename(UnitConverter) pism::units::Converter;
%shared_ptr(pism::units::System);
%feature("valuewrapper") pism::units::System;
%feature("valuewrapper") pism::units::Unit;
%include "util/Units.hh"

%shared_ptr(pism::MaxTimestep)
%include "util/MaxTimestep.hh"

%include pism_DM.i
%include pism_Vec.i
/* End of independent PISM classes. */

%shared_ptr(pism::Config);
%shared_ptr(pism::NetCDFConfig);
%shared_ptr(pism::DefaultConfig);
%include "util/ConfigInterface.hh"
%include "util/Config.hh"

#if (Pism_USE_JANSSON==1)
%shared_ptr(pism::ConfigJSON);
%include "util/ConfigJSON.hh"
#endif

/* EnthalpyConverter uses Config, so we need to wrap Config first (see above). */
%shared_ptr(pism::ColdEnthalpyConverter);
pism_class(pism::EnthalpyConverter, "pism/util/EnthalpyConverter.hh");

pism_class(pism::Time, "pism/util/Time.hh")

%include "util/Profiling.hh"
%shared_ptr(pism::Context);
%include "util/Context.hh"

%include pism_Grid.i

/* File uses Grid, so Grid has to be wrapped first. */
%include pism_File.i

/* make sure pism_File.i is included before VariableMetadata.hh */
%include pism_VariableMetadata.i

/* array::Array uses Grid and VariableMetadata so they have to be wrapped first. */
%include pism_Array.i

/* pism::Vars uses array::Array, so Array has to be wrapped first. */
%include pism_Vars.i


%shared_ptr(pism::Diagnostic)
%include "util/Diagnostic.hh"
%include "stressbalance/timestepping.hh"

%shared_ptr(pism::Component)
%include "util/Component.hh"

/* GeometryEvolution is a Component, so this has to go after Component.hh */
%include geometry.i

%include "basalstrength/basal_resistance.hh"

%include pism_FlowLaw.i

%include pism_ColumnSystem.i

%include pism_energy.i

/* SSAForwardRunFromInputFile sets up a yield stress model, which
 * requires a hydrology model.
 */
%include pism_Hydrology.i

%include "geometry/grounded_cell_fraction.hh"
%ignore pism::cell_type::UNKNOWN;
%include "util/cell_type.hh"
%include "pism_python.hh"

pism_class(pism::MohrCoulombPointwise, "pism/basalstrength/MohrCoulombPointwise.hh")
pism_class(pism::YieldStress, "pism/basalstrength/YieldStress.hh")
pism_class(pism::ConstantYieldStress, "pism/basalstrength/ConstantYieldStress.hh")
pism_class(pism::MohrCoulombYieldStress, "pism/basalstrength/MohrCoulombYieldStress.hh")
pism_class(pism::RegionalYieldStress, "pism/regional/RegionalYieldStress.hh")

%rename(StressBalanceInputs) pism::stressbalance::Inputs;

%include pism_SSA.i

%include pism_SIA.i

%include pism_blatter.i

%include pism_BedDef.i

%include AgeModel.i

/* The regional model implements some classes derived from SSAFD and
 * SIAFD, so this %include has to appear after %including the rest of
 * PISM's stress balance headers.
 */
%{
#include "regional/SSAFD_Regional.hh"
#include "regional/SIAFD_Regional.hh"
%}
%shared_ptr(pism::stressbalance::SSAFD_Regional)
%include "regional/SSAFD_Regional.hh"
%shared_ptr(pism::stressbalance::SIAFD_Regional)
%include "regional/SIAFD_Regional.hh"

%include "util/projection.hh"

%rename(linear_chi) pism::fem::linear::chi;
%rename(linear_n_chi) pism::fem::linear::n_chi;
%rename(q0_chi) pism::fem::q0::chi;
%rename(q0_n_chi) pism::fem::q0::n_chi;
%rename(q0_n_sides) pism::fem::q0::n_sides;
%rename(q1_chi) pism::fem::q1::chi;
%rename(q1_n_chi) pism::fem::q1::n_chi;
%rename(q1_n_sides) pism::fem::q1::n_sides;
%rename(Q1ElementGeometry) pism::fem::q1::ElementGeometry;
%rename(Q1BoundaryQuadrature2) pism::fem::q1::BoundaryQuadrature2;
%rename(p1_chi) pism::fem::p1::chi;
%rename(p1_n_chi) pism::fem::p1::n_chi;
%rename(p1_n_sides) pism::fem::p1::n_sides;
%rename(P1ElementGeometry) pism::fem::p1::ElementGeometry;
%rename(P1BoundaryQuadrature2) pism::fem::p1::BoundaryQuadrature2;
%rename(q13d_chi) pism::fem::q13d::chi;

%include "util/fem/FEM.hh"
%include "util/fem/DirichletData.hh"
%include "util/fem/Element.hh"
%include "util/fem/ElementIterator.hh"
%include "util/fem/Quadrature.hh"

%include "util/node_types.hh"

%include pism_inverse.i

%include "coupler/util/PCFactory.hh"

pism_class(pism::ForcingOptions, "pism/coupler/util/options.hh")

pism_class(pism::ScalarForcing, "pism/util/ScalarForcing.hh")

%shared_ptr(pism::PCFactory< pism::surface::SurfaceModel >)
%template(_SurfaceFactoryBase) pism::PCFactory<pism::surface::SurfaceModel>;

%shared_ptr(pism::PCFactory<pism::ocean::OceanModel>)
%template(_OceanFactoryBase) pism::PCFactory<pism::ocean::OceanModel>;

%shared_ptr(pism::PCFactory<pism::ocean::sea_level::SeaLevel>)
%template(_SeaLevelFactoryBase) pism::PCFactory<pism::ocean::sea_level::SeaLevel>;

%shared_ptr(pism::PCFactory<pism::atmosphere::AtmosphereModel>)
%template(_AtmosphereFactoryBase) pism::PCFactory<pism::atmosphere::AtmosphereModel>;

%include pism_ocean.i

%include pism_frontalmelt.i

/* surface models use atmosphere models as inputs so we need to define atmosphere models first */
%include pism_atmosphere.i

%include pism_surface.i

%include pism_calving.i

%include pism_verification.i

%include "energy/bootstrapping.hh"

#if (Pism_DEBUG==1)
pism_class(pism::Poisson, "pism/util/Poisson.hh")
#endif

pism_class(pism::FractureDensity, "pism/fracturedensity/FractureDensity.hh")
%include "util/label_components.hh"

pism_class(pism::IceModel, "pism/icemodel/IceModel.hh")
