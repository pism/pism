// Copyright (C) 2011, 2014, 2015, 2017, 2018 PISM Authors
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

#ifndef _PSFACTORY_H_
#define _PSFACTORY_H_

#include "pism/coupler/SurfaceModel.hh"
#include "pism/coupler/util/PCFactory.hh"

namespace pism {
namespace surface {
class Factory : public PCFactory<SurfaceModel> {
  typedef atmosphere::AtmosphereModel InputModel;

public:
  Factory(IceGrid::ConstPtr g, std::shared_ptr<InputModel> input);
  ~Factory();

  std::shared_ptr<SurfaceModel> create();
  std::shared_ptr<SurfaceModel> create(const std::string &type);

  void set_default(const std::string &type);
private:
  class SurfaceModelCreator {
  public:
    virtual std::shared_ptr<SurfaceModel> create(IceGrid::ConstPtr grid,
                                                 std::shared_ptr<InputModel> input) = 0;
    virtual ~SurfaceModelCreator() {}
  };

  // Creator for a specific model class M.
  template <class M>
  class SpecificSurfaceModelCreator : public SurfaceModelCreator {
  public:
    std::shared_ptr<SurfaceModel> create(IceGrid::ConstPtr grid,
                                         std::shared_ptr<InputModel> input) {
      return std::shared_ptr<SurfaceModel>(new M(grid, input));
    }
  };

  template <class M>
  void add_surface_model(const std::string &name) {
    m_surface_models[name].reset(new SpecificSurfaceModelCreator<M>);
  }

  std::shared_ptr<SurfaceModel> surface_model(const std::string &type,
                                              std::shared_ptr<InputModel> input);

  std::shared_ptr<InputModel> m_input;

  std::map<std::string, std::shared_ptr<SurfaceModelCreator> > m_surface_models;
};
} // end of namespace surface
} // end of namespace pism

#endif /* _PSFACTORY_H_ */
