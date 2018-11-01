%{
#include "coupler/ocean/Constant.hh"
#include "coupler/ocean/Cache.hh"
#include "coupler/ocean/ConstantPIK.hh"
#include "coupler/ocean/Delta_SMB.hh"
#include "coupler/ocean/Anomaly.hh"
#include "coupler/ocean/Delta_T.hh"
#include "coupler/ocean/Frac_SMB.hh"
#include "coupler/ocean/Frac_MBP.hh"
#include "coupler/ocean/GivenClimate.hh"
#include "coupler/ocean/GivenTH.hh"
#include "coupler/ocean/Pico.hh"
#include "coupler/ocean/PicoGeometry.hh"
#include "coupler/ocean/PicoPhysics.hh"
#include "coupler/ocean/sea_level/Delta_SL.hh"
#include "coupler/ocean/sea_level/Delta_SL_2D.hh"
%}

%shared_ptr(pism::ocean::OceanModel)
%include "coupler/OceanModel.hh"

%shared_ptr(pism::ocean::CompleteOceanModel)
%include "coupler/ocean/CompleteOceanModel.hh"

%shared_ptr(pism::ocean::Constant)
%rename(OceanConstant) pism::ocean::Constant;
%include "coupler/ocean/Constant.hh"

%shared_ptr(pism::ocean::Cache)
%rename(OceanCache) pism::ocean::Cache;
%include "coupler/ocean/Cache.hh"

%shared_ptr(pism::ocean::PIK)
%rename(OceanPIK) pism::ocean::PIK;
%include "coupler/ocean/ConstantPIK.hh"

%shared_ptr(pism::ocean::Delta_SMB)
%rename(OceanDeltaSMB) pism::ocean::Delta_SMB;
%include "coupler/ocean/Delta_SMB.hh"

%shared_ptr(pism::ocean::Anomaly)
%rename(OceanAnomaly) pism::ocean::Anomaly;
%include "coupler/ocean/Anomaly.hh"

%shared_ptr(pism::ocean::Delta_T)
%rename(OceanDeltaT) pism::ocean::Delta_T;
%include "coupler/ocean/Delta_T.hh"

%shared_ptr(pism::ocean::Frac_SMB)
%rename(OceanFracSMB) pism::ocean::Frac_SMB;
%include "coupler/ocean/Frac_SMB.hh"

%shared_ptr(pism::ocean::Frac_MBP)
%rename(OceanFracMBP) pism::ocean::Frac_MBP;
%include "coupler/ocean/Frac_MBP.hh"

%shared_ptr(pism::ocean::Given)
%rename(OceanGiven) pism::ocean::Given;
%include "coupler/ocean/GivenClimate.hh"

%shared_ptr(pism::ocean::GivenTH)
%rename(OceanGivenTH) pism::ocean::GivenTH;
%include "coupler/ocean/GivenTH.hh"

%shared_ptr(pism::ocean::Pico)
%rename(OceanPico) pism::ocean::Pico;
%include "coupler/ocean/Pico.hh"

%shared_ptr(pism::ocean::PicoGeometry)
%include "coupler/ocean/PicoGeometry.hh"

%shared_ptr(pism::ocean::PicoPhysics)
%include "coupler/ocean/PicoPhysics.hh"

%shared_ptr(pism::ocean::sea_level::SeaLevel)
%include "coupler/SeaLevel.hh"

%shared_ptr(pism::ocean::sea_level::Delta_SL)
%rename(SeaLevelDelta) pism::ocean::sea_level::Delta_SL;
%include "coupler/ocean/sea_level/Delta_SL.hh"

%shared_ptr(pism::ocean::sea_level::Delta_SL_2D)
%rename(SeaLevelDelta2D) pism::ocean::sea_level::Delta_SL_2D;
%include "coupler/ocean/sea_level/Delta_SL_2D.hh"
