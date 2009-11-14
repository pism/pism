// Copyright (C) 2008-2009 Ed Bueler, Constantine Khroulev, Gudfinna Adalgeirsdottir,
//                         and Andy Aschwanden
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef __pGreenlandAtmosCoupler_hh
#define __pGreenlandAtmosCoupler_hh

#include <petsc.h>
#include "../base/iceModelVec.hh"
#include "../base/PISMVars.hh"
#include "localMassBalance.hh"
#include "iceModelVec2T.hh"
#include "pccoupler.hh"


//! A derived class of PISMAtmosphereCoupler which has a simple snow process model suitable for Greenland modeling.
/*!
There are two parts to a "snow process model".
-# either a (snow/surface) temperature parameterization or stored (snow/surface)
   temperature maps,
-# a choice of PDD models, which convert the number of positive degree days, and
   the snow precipitation rate, into a surface mass balance (flux).

The default temperature parameterization is from \ref Faustoetal2009, directly
applicable to Greenland.  The constants can be changed, however, to describe other 
parameterizations which depend linearly on the same things (latitude, longitude,
and surface elevation).

If stored snow-surface temperature maps are used they are read through 
IceModelVec2T *snowtemps.

The PDD schemes are accessed through LocalMassBalance *mbscheme.
The default PDD model uses formula (6) from \ref Faustoetal2009, and uses a 
\ref CalovGreve05 evaluation of the expectation integral for the number of positive
degree days, given a normally-distributed daily variability.
An alternative PDD scheme simulates the daily variability directly by a normal
(pseudo-) random variable.

This model makes a severe simplifying assumption about temperatures.  Namely, 
vsurftemp is \e both the ice upper boundary temperature (i.e. below completion
of firn processes ...) \e and it is used as the mean annual temperature 
used in the snow processes model.  That is, the temperature used in the 
snow processes model to convert the snow precipitation rate vsnowprecip into 
surface mass flux vsurfmassflux is here assumed to be equal to vsurftemp.
In fact, the available temperature in real modeling is frequently the +2 m
temperature, in the air above the snow surface at the level measured by a met
station, and/or if may be from an atmospheric circulation model.  Thus there
are three temperatures we might care about, all roughly at the surface of the ice
sheet:

-# the available 2m temperature (air)
-# the needed temperature for correctly modeling melting in the mass balance scheme, and
-# the needed temperature for the upper boundary of the conservation of energy scheme.

Here we just assume they are all the same.  It is a deficiency of the model.
 */
class PISMGreenlandAtmosCoupler : public PISMAtmosphereCoupler {

public:
  PISMGreenlandAtmosCoupler();
  virtual ~PISMGreenlandAtmosCoupler();

  using PISMAtmosphereCoupler::initFromOptions;
  virtual PetscErrorCode initFromOptions(IceGrid* g, const PISMVars &variables);

  virtual PetscErrorCode setLMBScheme(LocalMassBalance *usethisscheme);

  using PISMAtmosphereCoupler::writeCouplingFieldsToFile;
  virtual PetscErrorCode writeCouplingFieldsToFile(
             PetscScalar t_years, const char *filename);

  using PISMAtmosphereCoupler::updateSurfMassFluxAndProvide;
  virtual PetscErrorCode updateSurfMassFluxAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvsmf);

  using PISMAtmosphereCoupler::updateSurfTempAndProvide;
  virtual PetscErrorCode updateSurfTempAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvst);

  virtual PetscErrorCode max_timestep(PetscScalar t_years, PetscScalar &dt_years);
protected:
  //! Defaults to the \ref Faustoetal2009 scheme.  Called when no stored temperature maps are available.
  /*!
    Computes the mean annual temperature as function of latitude, longitude, 
    surface elevation, and so on.
    Depends on the current state of IceModel fields, through pointer info.
   */
  virtual PetscErrorCode parameterizedUpdateSnowSurfaceTemp(
              PetscScalar t_years, PetscScalar dt_years);

  //! Instead of temperature parameterization we can use stored temperature maps read from a file.
  IceModelVec2T *snowtempmaps;

  LocalMassBalance *mbscheme;

  //! The snow precipitation rate in ice-equivalent meters per second.
  /*! vsurfmassflux is computed by LocalMassBalance scheme from this rate.  */
  IceModelVec2 vsnowprecip;

  //! The snow precipitation rate in ice-equivalent meters per second, read from a file.
  IceModelVec2T *snowprecipmaps;

  //! The mean July (julian day = 196; July 15) snow temperature used in the mass balance scheme.
  /*! This field is a diagnostic extra output of PISMGreenlandAtmosCoupler.  IceModel never
      gets a pointer to it. */
  IceModelVec2 vsnowtemp_mj;  

  //! Pointers to IceModelVecs provided by IceModel.
  IceModelVec2 *surfelev, *lat, *lon;
};

#endif

