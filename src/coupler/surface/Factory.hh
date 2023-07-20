// Copyright (C) 2011, 2014, 2015, 2017, 2018, 2019, 2021, 2023 PISM Authors
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

#ifndef PISM_SURFACE_FACTORY
#define PISM_SURFACE_FACTORY

#include "pism/coupler/SurfaceModel.hh"
#include "pism/coupler/util/PCFactory.hh"

namespace pism {
namespace surface {

class Factory : public PCFactory<SurfaceModel> {
public:
  typedef std::shared_ptr<atmosphere::AtmosphereModel> AtmospherePtr;

  Factory(std::shared_ptr<const Grid> g, AtmospherePtr input);

  ~Factory() = default;

private:
  // Creator for a specific model class M.
  template <class M>
  class SurfaceModelCreator : public ModelCreator {
  public:
    SurfaceModelCreator(AtmospherePtr input) : m_input(input) {
      // empty
    }

    std::shared_ptr<SurfaceModel> create(std::shared_ptr<const Grid> grid) {
      return std::make_shared<M>(grid, m_input);
    }

  private:
    AtmospherePtr m_input;
  };

  template <class M>
  void add_surface_model(const std::string &name) {
    PCFactory<SurfaceModel>::m_models[name] = std::make_shared<SurfaceModelCreator<M> >(m_input);
  }

  AtmospherePtr m_input;
};
} // end of namespace surface
} // end of namespace pism

#endif /* PISM_SURFACE_FACTORY */
