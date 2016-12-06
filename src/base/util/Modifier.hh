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
#include "base/util/pism_utilities.hh"

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
    : Model(g), m_input_model(in) {}
  virtual ~Modifier()
  {
    if (m_input_model != NULL) {
      delete m_input_model;
    }
  }

protected:
  virtual MaxTimestep max_timestep_impl(double my_t) const
  {
    if (m_input_model != NULL) {
      return m_input_model->max_timestep(my_t);
    } else {
      return MaxTimestep();
    }
  }

  virtual void update_impl(double t, double dt)
  {
    Model::m_t = t;
    Model::m_dt = dt;
    if (m_input_model != NULL) {
      m_input_model->update(t, dt);
    }
  }

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const
  {
    // give the model class a chance to add diagnostics
    std::map<std::string, Diagnostic::Ptr> result = Model::diagnostics_impl();

    // add diagnostics from an input model, if it exists
    if (m_input_model != NULL) {
      result = pism::combine(result, m_input_model->diagnostics());
    }
    return result;
  }

  virtual std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const
  {
    // give the model class a chance to add diagnostics
    std::map<std::string, TSDiagnostic::Ptr> result = Model::ts_diagnostics_impl();

    // add diagnostics from an input model, if it exists
    if (m_input_model != NULL) {
      result = pism::combine(result, m_input_model->ts_diagnostics());
    }
    return result;
  }

  virtual void define_model_state_impl(const PIO &output) const {
    if (m_input_model != NULL) {
      m_input_model->define_model_state(output);
    }
  }

  virtual void write_model_state_impl(const PIO &output) const {
    if (m_input_model != NULL) {
      m_input_model->write_model_state(output);
    }
  }
protected:
  Model *m_input_model;
};  
} // end of namespace pism

#endif /* _MODIFIER_H_ */
