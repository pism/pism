// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#ifndef _PISMMOHRCOULOMBYIELDSTRESS_H_
#define _PISMMOHRCOULOMBYIELDSTRESS_H_

#include "PISMDiagnostic.hh"
#include "PISMYieldStress.hh"
#include "PISMHydrology.hh"
#include "iceModelVec.hh"

namespace pism {


//! \brief PISM's default basal yield stress model which applies the Mohr-Coulomb model of deformable, pressurized till.
class MohrCoulombYieldStress : public YieldStress
{
public:
  MohrCoulombYieldStress(IceGrid &g, const Config &conf, Hydrology *hydro);

  virtual ~MohrCoulombYieldStress();

  virtual void init(Vars &vars);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);

  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);

  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);

  virtual void update(double my_t, double my_dt);

  virtual void basal_material_yield_stress(IceModelVec2S &result);

protected:
  IceModelVec2S m_till_phi, m_tauc, m_tillwat, m_Po;
  IceModelVec2S m_bwat;  // only allocated and used if tauc_add_transportable_water = true
  IceModelVec2S *m_bed_topography;
  IceModelVec2Int *m_mask;
  Vars *m_variables;
  Hydrology *m_hydrology;

  void topg_to_phi();
  void tauc_to_phi();
};

} // end of namespace pism

#endif /* _PISMMOHRCOULOMBYIELDSTRESS_H_ */
