/* Copyright (C) 2015, 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _MODIFIER_H_
#define _MODIFIER_H_

#include "base/util/MaxTimestep.hh"

namespace pism {
//! \brief This template allows creating Component_TS (AtmosphereModel,
//! SurfaceModel and OceanModel) modifiers with minimum effort.
/*!
 * A specialization of this template will implement all important methods
 * except init(). This means that to create a complete modifier, one needs to
 * re-implement interesting methods, without worrying about preserving
 * modifier's "transparency".
 */
template<class Model>
class Modifier : public Model
{
public:
  Modifier(IceGrid::ConstPtr g, Model* in)
    : Model(g), input_model(in) {}
  virtual ~Modifier()
  {
    if (input_model != NULL) {
      delete input_model;
    }
  }

protected:
  virtual MaxTimestep max_timestep_impl(double my_t)
  {
    if (input_model != NULL) {
      return input_model->max_timestep(my_t);
    } else {
      return MaxTimestep();
    }
  }

  virtual void update_impl(double my_t, double my_dt)
  {
    Model::m_t = my_t;
    Model::m_dt = my_dt;
    if (input_model != NULL) {
      input_model->update(my_t, my_dt);
    }
  }

  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {
    if (input_model != NULL) {
      input_model->get_diagnostics(dict, ts_dict);
    }
  }

  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc)
  {
    if (input_model != NULL) {
      input_model->write_variables(vars, nc);
    }
  }

  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result)
  {
    if (input_model != NULL) {
      input_model->add_vars_to_output(keyword, result);
    }
  }

  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                     IO_Type nctype) {
    if (input_model != NULL) {
      input_model->define_variables(vars, nc, nctype);
    }
  }
protected:
  Model *input_model;
};  
} // end of namespace pism

#endif /* _MODIFIER_H_ */
