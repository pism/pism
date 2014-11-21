// Copyright (C) 2012-2014 PISM Authors
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


#include "hydrology_diagnostics.hh"

namespace pism {

Hydrology_bwat::Hydrology_bwat(Hydrology *m, IceGrid &g, Vars &my_vars)
    : Diag<Hydrology>(m, g, my_vars) {
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "bwat", grid));
  set_attrs("thickness of transportable water in subglacial layer", "", "m", "m", 0);
}

void Hydrology_bwat::compute(IceModelVec* &output) {
  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "bwat", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  model->subglacial_water_thickness(*result);
  output = result;
}

Hydrology_bwp::Hydrology_bwp(Hydrology *m, IceGrid &g, Vars &my_vars)
    : Diag<Hydrology>(m, g, my_vars) {
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "bwp", grid));
  set_attrs("pressure of transportable water in subglacial layer", "", "Pa", "Pa", 0);
}


void Hydrology_bwp::compute(IceModelVec* &output) {
  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "bwp", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  model->subglacial_water_pressure(*result);
  output = result;
}


Hydrology_bwprel::Hydrology_bwprel(Hydrology *m, IceGrid &g, Vars &my_vars)
    : Diag<Hydrology>(m, g, my_vars) {
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "bwprel", grid));
  set_attrs("pressure of transportable water in subglacial layer as fraction of the overburden pressure", "",
            "", "", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}


void Hydrology_bwprel::compute(IceModelVec* &output) {
  double fill = grid.config.get("fill_value");
  IceModelVec2S *Po     = new IceModelVec2S,
                *result = new IceModelVec2S;
  result->create(grid, "bwprel", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  Po->create(grid, "Po_temporary", WITHOUT_GHOSTS);

  model->subglacial_water_pressure(*result);
  model->overburden_pressure(*Po);

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*Po);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*Po)(i,j) > 0.0) {
      (*result)(i,j) /= (*Po)(i,j);
    } else {
      (*result)(i,j) = fill;
    }
  }

  output = result;
}


Hydrology_effbwp::Hydrology_effbwp(Hydrology *m, IceGrid &g, Vars &my_vars)
    : Diag<Hydrology>(m, g, my_vars) {
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "effbwp", grid));
  set_attrs("effective pressure of transportable water in subglacial layer (overburden pressure minus water pressure)",
            "", "Pa", "Pa", 0);
}


void Hydrology_effbwp::compute(IceModelVec* &output) {
  IceModelVec2S *P      = new IceModelVec2S,
                *result = new IceModelVec2S;
  result->create(grid, "effbwp", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  P->create(grid, "P_temporary", WITHOUT_GHOSTS);

  model->subglacial_water_pressure(*P);
  model->overburden_pressure(*result);
  result->add(-1.0,*P);  // result <-- result + (-1.0) P = Po - P

  output = result;
}


Hydrology_hydrobmelt::Hydrology_hydrobmelt(Hydrology *m, IceGrid &g, Vars &my_vars)
    : Diag<Hydrology>(m, g, my_vars) {
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "hydrobmelt", grid));
  set_attrs("the version of bmelt seen by the hydrology model",
            "", "m s-1", "m/year", 0);
}


void Hydrology_hydrobmelt::compute(IceModelVec* &output) {
  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "hydrobmelt", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  // the value reported diagnostically is merely the last value filled
  (model->bmelt_local).copy_to(*result);
  output = result;
}


Hydrology_hydroinput::Hydrology_hydroinput(Hydrology *m, IceGrid &g, Vars &my_vars)
    : Diag<Hydrology>(m, g, my_vars) {
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "hydroinput", grid));
  set_attrs("total water input into subglacial hydrology layer",
            "", "m s-1", "m/year", 0);
}


void Hydrology_hydroinput::compute(IceModelVec* &output) {
  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "hydroinput", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  // the value reported diagnostically is merely the last value filled
  (model->total_input).copy_to(*result);
  output = result;
}


Hydrology_wallmelt::Hydrology_wallmelt(Hydrology *m, IceGrid &g, Vars &my_vars)
    : Diag<Hydrology>(m, g, my_vars) {
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "wallmelt", grid));
  set_attrs("wall melt into subglacial hydrology layer from (turbulent) dissipation of energy in transportable water",
            "", "m s-1", "m/year", 0);
}


void Hydrology_wallmelt::compute(IceModelVec* &output) {
  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "wallmelt", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  model->wall_melt(*result);
  output = result;
}


} // end of namespace pism
