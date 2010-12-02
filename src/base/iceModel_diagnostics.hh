// Copyright (C) 2010 Constantine Khroulev
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

#ifndef _ICEMODEL_DIAGNOSTICS_H_
#define _ICEMODEL_DIAGNOSTICS_H_

#include "iceModel.hh"

//! \brief Computes vertically-averaged ice hardness.
class IceModel_hardav : public PISMDiag<IceModel>
{
public:
  IceModel_hardav(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes a diagnostic field filled with processor rank values.
class IceModel_rank : public PISMDiag<IceModel>
{
public:
  IceModel_rank(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes basal (pore) water pressure.
class IceModel_bwp : public PISMDiag<IceModel>
{
public:
  IceModel_bwp(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes CTS, CTS = E/E_s(p).
class IceModel_cts : public PISMDiag<IceModel>
{
public:
  IceModel_cts(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the rate of change of ice surface elevation as a sum of the
//! bedrock uplift rate and the thickness rate of change.
class IceModel_dhdt : public PISMDiag<IceModel>
{
public:
  IceModel_dhdt(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the number of ice-filled cells is a processor's domain.
class IceModel_proc_ice_area : public PISMDiag<IceModel>
{
public:
  IceModel_proc_ice_area(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes ice temperature from enthalpy.
class IceModel_temp : public PISMDiag<IceModel>
{
public:
  IceModel_temp(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Compute the pressure-adjusted temperature in degrees C corresponding
//! to ice temperature.
class IceModel_temp_pa : public PISMDiag<IceModel>
{
public:
  IceModel_temp_pa(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes basal values of the pressure-adjusted temperature.
class IceModel_temppabase : public PISMDiag<IceModel>
{
public:
  IceModel_temppabase(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes surface values of ice enthalpy.
class IceModel_enthalpysurf : public PISMDiag<IceModel>
{
public:
  IceModel_enthalpysurf(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes enthalpy at the base of the ice.
class IceModel_enthalpybase : public PISMDiag<IceModel>
{
public:
  IceModel_enthalpybase(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes ice temperature at the base of the ice.
class IceModel_tempbase : public PISMDiag<IceModel>
{
public:
  IceModel_tempbase(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes ice temperature at the surface of the ice.
class IceModel_tempsurf : public PISMDiag<IceModel>
{
public:
  IceModel_tempsurf(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the liquid water fraction.
class IceModel_liqfrac : public PISMDiag<IceModel>
{
public:
  IceModel_liqfrac(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the total thickness of temperate ice in a column.
class IceModel_tempicethk : public PISMDiag<IceModel>
{
public:
  IceModel_tempicethk(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};
//! \brief Computes the thickness of the basal layer of temperate ice.
class IceModel_tempicethk_basal : public PISMDiag<IceModel>
{
public:
  IceModel_tempicethk_basal(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

#endif /* _ICEMODEL_DIAGNOSTICS_H_ */

