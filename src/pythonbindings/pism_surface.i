%{
#include "coupler/surface/EISMINTII.hh"
#include "coupler/surface/Delta_T.hh"
#include "coupler/surface/ConstantPIK.hh"
#include "coupler/surface/Cache.hh"
#include "coupler/surface/Anomaly.hh"
#include "coupler/surface/Elevation.hh"
#include "coupler/surface/LapseRates.hh"
#include "coupler/surface/Simple.hh"
#include "coupler/surface/TemperatureIndex.hh"
#include "coupler/surface/GivenClimate.hh"
%}

%shared_ptr(pism::surface::SurfaceModel)
%include "coupler/SurfaceModel.hh"

%shared_ptr(pism::surface::Delta_T)
%rename(SurfaceDeltaT) pism::surface::Delta_T;
%include "coupler/surface/Delta_T.hh"

%shared_ptr(pism::surface::PSFormulas)
%include "coupler/surface/Formulas.hh"

%shared_ptr(pism::surface::EISMINTII)
%rename(SurfaceEISMINTII) pism::surface::EISMINTII;
%feature("notabstract") pism::surface::EISMINTII;
%include "coupler/surface/EISMINTII.hh"

%shared_ptr(pism::surface::PIK)
%rename(SurfacePIK) pism::surface::PIK;
%include "coupler/surface/ConstantPIK.hh"

%shared_ptr(pism::surface::Cache)
%rename(SurfaceCache) pism::surface::Cache;
%include "coupler/surface/Cache.hh"

%shared_ptr(pism::surface::Anomaly)
%rename(SurfaceAnomaly) pism::surface::Anomaly;
%include "coupler/surface/Anomaly.hh"

%shared_ptr(pism::surface::Given)
%rename(SurfaceGiven) pism::surface::Given;
%include "coupler/surface/GivenClimate.hh"

%shared_ptr(pism::surface::Elevation)
%rename(SurfaceElevation) pism::surface::Elevation;
%include "coupler/surface/Elevation.hh"

%shared_ptr(pism::surface::LapseRates)
%rename(SurfaceLapseRates) pism::surface::LapseRates;
%include "coupler/surface/LapseRates.hh"

%shared_ptr(pism::surface::Simple)
%rename(SurfaceSimple) pism::surface::Simple;
%include "coupler/surface/Simple.hh"

%shared_ptr(pism::surface::TemperatureIndex)
%rename(SurfaceTemperatureIndex) pism::surface::TemperatureIndex;
%include "coupler/surface/TemperatureIndex.hh"
