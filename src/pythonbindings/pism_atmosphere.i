%{
#include "coupler/atmosphere/Anomaly.hh"
#include "coupler/atmosphere/CosineYearlyCycle.hh"
#include "coupler/atmosphere/Delta_P.hh"
#include "coupler/atmosphere/Delta_T.hh"
#include "coupler/atmosphere/Factory.hh"
#include "coupler/atmosphere/Frac_P.hh"
#include "coupler/atmosphere/GivenClimate.hh"
#include "coupler/atmosphere/ElevationChange.hh"
#include "coupler/atmosphere/PIK.hh"
#include "coupler/atmosphere/PrecipitationScaling.hh"
#include "coupler/atmosphere/SeariseGreenland.hh"
#include "coupler/atmosphere/Uniform.hh"
#include "coupler/atmosphere/WeatherStation.hh"
#include "coupler/atmosphere/YearlyCycle.hh"
#include "coupler/atmosphere/OrographicPrecipitation.hh"
%}

%shared_ptr(pism::atmosphere::AtmosphereModel)
%include "coupler/AtmosphereModel.hh"

%shared_ptr(pism::atmosphere::Anomaly)
%rename(AtmosphereAnomaly) pism::atmosphere::Anomaly;
%include "coupler/atmosphere/Anomaly.hh"

%shared_ptr(pism::atmosphere::Delta_P)
%rename(AtmosphereDeltaP) pism::atmosphere::Delta_P;
%include "coupler/atmosphere/Delta_P.hh"

%shared_ptr(pism::atmosphere::Delta_T)
%rename(AtmosphereDeltaT) pism::atmosphere::Delta_T;
%include "coupler/atmosphere/Delta_T.hh"

%shared_ptr(pism::atmosphere::Frac_P)
%rename(AtmosphereFracP) pism::atmosphere::Frac_P;
%include "coupler/atmosphere/Frac_P.hh"

%shared_ptr(pism::atmosphere::Given)
%rename(AtmosphereGiven) pism::atmosphere::Given;
%include "coupler/atmosphere/GivenClimate.hh"

%shared_ptr(pism::atmosphere::ElevationChange)
%rename(AtmosphereElevationChange) pism::atmosphere::ElevationChange;
%include "coupler/atmosphere/ElevationChange.hh"

%shared_ptr(pism::atmosphere::PrecipitationScaling)
%rename(AtmospherePrecipScaling) pism::atmosphere::PrecipitationScaling;
%include "coupler/atmosphere/PrecipitationScaling.hh"

%shared_ptr(pism::atmosphere::WeatherStation)
%rename(AtmosphereWeatherStation) pism::atmosphere::WeatherStation;
%include "coupler/atmosphere/WeatherStation.hh"

%shared_ptr(pism::atmosphere::YearlyCycle)
%rename(AtmosphereYearlyCycle) pism::atmosphere::YearlyCycle;
%include "coupler/atmosphere/YearlyCycle.hh"

%shared_ptr(pism::atmosphere::CosineYearlyCycle)
%rename(AtmosphereCosineYearlyCycle) pism::atmosphere::CosineYearlyCycle;
%include "coupler/atmosphere/CosineYearlyCycle.hh"

%shared_ptr(pism::atmosphere::PIK)
%rename(AtmospherePIK) pism::atmosphere::PIK;
%include "coupler/atmosphere/PIK.hh"

%shared_ptr(pism::atmosphere::SeaRISEGreenland)
%rename(AtmosphereSeaRISEGreenland) pism::atmosphere::SeaRISEGreenland;
%include "coupler/atmosphere/SeariseGreenland.hh"

%shared_ptr(pism::atmosphere::Uniform)
%rename(AtmosphereUniform) pism::atmosphere::Uniform;
%include "coupler/atmosphere/Uniform.hh"

%shared_ptr(pism::atmosphere::Factory)
%rename(AtmosphereFactory) pism::atmosphere::Factory;
%include "coupler/atmosphere/Factory.hh"

%shared_ptr(pism::atmosphere::OrographicPrecipitation)
%rename(AtmosphereOrographicPrecipitation) pism::atmosphere::OrographicPrecipitation;
%include "coupler/atmosphere/OrographicPrecipitation.hh"
