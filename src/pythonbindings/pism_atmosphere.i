%{
#include "coupler/atmosphere/Anomaly.hh"
#include "coupler/atmosphere/ConstantPIK.hh"
#include "coupler/atmosphere/CosineYearlyCycle.hh"
#include "coupler/atmosphere/Delta_P.hh"
#include "coupler/atmosphere/Delta_T.hh"
#include "coupler/atmosphere/Factory.hh"
#include "coupler/atmosphere/Frac_P.hh"
#include "coupler/atmosphere/GivenClimate.hh"
#include "coupler/atmosphere/Paleo_precip.hh"
#include "coupler/atmosphere/SeariseGreenland.hh"
#include "coupler/atmosphere/Uniform.hh"
#include "coupler/atmosphere/WeatherStation.hh"
#include "coupler/atmosphere/YearlyCycle.hh"
%}

%shared_ptr(pism::atmosphere::AtmosphereModel)
%include "coupler/AtmosphereModel.hh"

%shared_ptr(pism::atmosphere::Anomaly)
%rename(AtmosphereAnomaly) pism::atmosphere::Anomaly;
%include "coupler/atmosphere/Anomaly.hh"

%shared_ptr(pism::atmosphere::PIK)
%rename(AtmospherePIK) pism::atmosphere::PIK;
%include "coupler/atmosphere/ConstantPIK.hh"

%shared_ptr(pism::atmosphere::Delta_P)
%rename(AtmosphereDelta_P) pism::atmosphere::Delta_P;
%include "coupler/atmosphere/Delta_P.hh"

%shared_ptr(pism::atmosphere::Delta_T)
%rename(AtmosphereDelta_T) pism::atmosphere::Delta_T;
%include "coupler/atmosphere/Delta_T.hh"

%shared_ptr(pism::atmosphere::Frac_P)
%rename(AtmosphereFrac_P) pism::atmosphere::Frac_P;
%include "coupler/atmosphere/Frac_P.hh"

%shared_ptr(pism::atmosphere::Given)
%rename(AtmosphereGiven) pism::atmosphere::Given;
%include "coupler/atmosphere/GivenClimate.hh"

%shared_ptr(pism::atmosphere::PaleoPrecip)
%rename(AtmospherePaleo_precip) pism::atmosphere::PaleoPrecip;
%include "coupler/atmosphere/Paleo_precip.hh"

%shared_ptr(pism::atmosphere::WeatherStation)
%rename(AtmosphereWeatherStation) pism::atmosphere::WeatherStation;
%include "coupler/atmosphere/WeatherStation.hh"

%shared_ptr(pism::atmosphere::YearlyCycle)
%rename(AtmosphereYearlyCycle) pism::atmosphere::YearlyCycle;
%include "coupler/atmosphere/YearlyCycle.hh"

%shared_ptr(pism::atmosphere::CosineYearlyCycle)
%rename(AtmosphereCosineYearlyCycle) pism::atmosphere::CosineYearlyCycle;
%include "coupler/atmosphere/CosineYearlyCycle.hh"

%shared_ptr(pism::atmosphere::SeaRISEGreenland)
%rename(AtmosphereSeaRISEGreenland) pism::atmosphere::SeaRISEGreenland;
%include "coupler/atmosphere/SeariseGreenland.hh"

%shared_ptr(pism::atmosphere::Uniform)
%rename(AtmosphereUniform) pism::atmosphere::Uniform;
%include "coupler/atmosphere/Uniform.hh"
