/* Copyright (C) 2015, 2016, 2017, 2018 PISM Authors
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

#include <memory>               // std::unique_ptr

#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
//! \brief This template allows creating Component (AtmosphereModel,
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
  virtual ~Modifier() {}

protected:
  virtual MaxTimestep max_timestep_impl(double my_t) const
  {
    if (m_input_model) {
      return m_input_model->max_timestep(my_t);
    } else {
      return MaxTimestep();
    }
  }

  virtual void update_impl(double t, double dt)
  {
    Model::m_t = t;
    Model::m_dt = dt;
    if (m_input_model) {
      m_input_model->update(t, dt);
    }
  }

  virtual DiagnosticList diagnostics_impl() const
  {
    // give the model class a chance to add diagnostics
    DiagnosticList result = Model::diagnostics_impl();

    // add diagnostics from an input model, if it exists
    if (m_input_model) {
      result = pism::combine(result, m_input_model->diagnostics());
    }
    return result;
  }

  virtual TSDiagnosticList ts_diagnostics_impl() const
  {
    // give the model class a chance to add diagnostics
    TSDiagnosticList result = Model::ts_diagnostics_impl();

    // add diagnostics from an input model, if it exists
    if (m_input_model) {
      result = pism::combine(result, m_input_model->ts_diagnostics());
    }
    return result;
  }

  virtual void define_model_state_impl(const PIO &output) const {
    if (m_input_model) {
      m_input_model->define_model_state(output);
    }
  }

  virtual void write_model_state_impl(const PIO &output) const {
    if (m_input_model) {
      m_input_model->write_model_state(output);
    }
  }
protected:
  std::unique_ptr<Model> m_input_model;
};  
} // end of namespace pism

#endif /* _MODIFIER_H_ */
