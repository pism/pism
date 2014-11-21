// Copyright (C) 2010, 2011, 2012, 2013, 2014 PISM Authors
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

#ifndef __BedDef_hh
#define __BedDef_hh

#include "PISMComponent.hh"
#include "iceModelVec.hh"

namespace pism {

//! PISM bed deformation model (base class).
/*! Unlike other Component_TS derived classes, the update() method of
  BedDef has side-effects (modifies IceModel data members).
*/
class BedDef : public Component_TS {
public:
  BedDef(IceGrid &g, const Config &conf);
  virtual ~BedDef() {}
  virtual void init(Vars &vars);
  virtual void update(double my_t, double my_dt) = 0;
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
protected:
  PetscErrorCode pismbeddef_allocate(); // packaged to simplify error checking
  void compute_uplift(double dt_beddef);
  double t_beddef_last;         //!< last bed deformation update year

  IceModelVec2S topg_initial;
  IceModelVec2S topg_last;      //!< last bed elevation
  IceModelVec2S *thk,           //!< pointer to the current ice thickness
    *topg,                      //!< pointer to the current bed elevation
    *uplift;                    //!< pointer to the bed uplift rate field
};

//! Pointwide isostasy bed deformation model.
class PBPointwiseIsostasy : public BedDef {
public:
  PBPointwiseIsostasy(IceGrid &g, const Config &conf); 
  virtual ~PBPointwiseIsostasy() {}
  virtual void init(Vars &vars);
  virtual void update(double my_t, double my_dt);
protected:
  PetscErrorCode allocate();
  IceModelVec2S thk_last;       //!< last ice thickness
};

} // end of namespace pism

#endif  // __BedDef_hh
