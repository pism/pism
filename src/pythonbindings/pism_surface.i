%{
#include "coupler/surface/ISMIP6Climate.hh"
#include "coupler/surface/EISMINTII.hh"
#include "coupler/surface/Delta_T.hh"
#include "coupler/surface/ConstantPIK.hh"
#include "coupler/surface/Cache.hh"
#include "coupler/surface/Anomaly.hh"
#include "coupler/surface/Elevation.hh"
#include "coupler/surface/ElevationChange.hh"
#include "coupler/surface/Simple.hh"
#include "coupler/surface/TemperatureIndex.hh"
#include "coupler/surface/GivenClimate.hh"
#include "coupler/surface/ForceThickness.hh"
#include "coupler/surface/Initialization.hh"
#include "coupler/surface/Factory.hh"
#include "coupler/surface/DEBMSimplePointwise.hh"
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

%shared_ptr(pism::surface::ISMIP6)
%rename(SurfaceISMIP6) pism::surface::ISMIP6;
%include "coupler/surface/ISMIP6Climate.hh"

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

%shared_ptr(pism::surface::ElevationChange)
%rename(SurfaceElevationChange) pism::surface::ElevationChange;
%include "coupler/surface/ElevationChange.hh"

%shared_ptr(pism::surface::Simple)
%rename(SurfaceSimple) pism::surface::Simple;
%include "coupler/surface/Simple.hh"

%shared_ptr(pism::surface::TemperatureIndex)
%rename(SurfaceTemperatureIndex) pism::surface::TemperatureIndex;
%include "coupler/surface/TemperatureIndex.hh"

%shared_ptr(pism::surface::ForceThickness)
%rename(SurfaceForceThickness) pism::surface::ForceThickness;
%include "coupler/surface/ForceThickness.hh"

%shared_ptr(pism::surface::InitializationHelper)
%rename(SurfaceInitialization) pism::surface::InitializationHelper;
%include "coupler/surface/Initialization.hh"

%shared_ptr(pism::surface::Factory)

%rename(SurfaceFactory) pism::surface::Factory;
%include "coupler/surface/Factory.hh"

%shared_ptr(pism::surface::DEBMSimplePointwise)
%rename(SurfaceDEBMSimplePointwise) pism::surface::DEBMSimplePointwise;
%include "coupler/surface/DEBMSimplePointwise.hh"
