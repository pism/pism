// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 PISM Authors
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
  BedDef(const IceGrid &g);
  virtual ~BedDef() {}
  virtual void init();
  virtual void update(double my_t, double my_dt) = 0;

  const IceModelVec2S& bed_elevation() const;
  const IceModelVec2S& uplift() const;

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
protected:
  PetscErrorCode pismbeddef_allocate(); // packaged to simplify error checking
  void compute_uplift(double dt_beddef);

  //! time of the last bed deformation update
  double m_t_beddef_last;

  //! current bed elevation
  IceModelVec2S m_topg;

  //! bed elevation at the beginning of a run
  IceModelVec2S m_topg_initial;

  //! bed elevation at the time of the last update
  IceModelVec2S m_topg_last;

  //! bed uplift rate
  IceModelVec2S m_uplift;

  //! pointer to the current ice thickness
  IceModelVec2S *m_thk;
};

//! Pointwide isostasy bed deformation model.
class PBPointwiseIsostasy : public BedDef {
public:
  PBPointwiseIsostasy(const IceGrid &g); 
  virtual ~PBPointwiseIsostasy();
  virtual void init();
  virtual void update(double my_t, double my_dt);
protected:
  PetscErrorCode allocate();
  IceModelVec2S m_thk_last;       //!< last ice thickness
};

} // end of namespace pism

#endif  // __BedDef_hh
