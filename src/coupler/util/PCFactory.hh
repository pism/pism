// Copyright (C) 2011, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _PCFACTORY_H_
#define _PCFACTORY_H_

#include <memory>
#include <map>

#include "pism/util/IceGrid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_options.hh"

namespace pism {

class Config;

template <class Model>
class PCFactory {
public:

  PCFactory<Model>(IceGrid::ConstPtr g)
  : m_grid(g) {}
  ~PCFactory<Model>() {}

  //! Creates a boundary model. Processes command-line options.
  virtual std::shared_ptr<Model> create() {
    // build a list of available models:
    auto model_list = key_list(m_models);

    // build a list of available modifiers:
    auto modifier_list = key_list(m_modifiers);

    std::string description = ("Sets up the PISM " + m_option + " model."
                               " Available models: " + model_list +
                               " Available modifiers: " + modifier_list);

    // Get the command-line option:
    options::StringList choices("-" + m_option, description, m_default_type);

    return create(choices.to_string());
  }

  //! Creates a boundary model.
  virtual std::shared_ptr<Model> create(const std::string &type) {

    std::vector<std::string> choices = split(type, ',');

    // the first element has to be an *actual* model (not a modifier)
    auto j = choices.begin();

    auto result = model(*j);

    ++j;

    // process remaining arguments:
    for (;j != choices.end(); ++j) {
      result = modifier(*j, result);
    }

    return result;
  }

  //! Sets the default type name.
  virtual void set_default(const std::string &name) {
    if (m_models.find(name) == m_models.end()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "type %s is not registered", name.c_str());
    } else {
      m_default_type = name;
    }
  }

protected:

  //! Adds a boundary model to the dictionary.
  template <class M>
  void add_model(const std::string &name) {
    m_models[name].reset(new SpecificModelCreator<M>);
  }

  template <class M>
  void add_modifier(const std::string &name) {
    m_modifiers[name].reset(new SpecificModifierCreator<M>);
  }

  template<typename T>
  std::string key_list(std::map<std::string, T> list) {
    std::string result;
    auto k = list.begin();
    result = "[" + (k++)->first;
    for (; k != list.end(); k++) {
      result += ", " + k->first;
    }
    result += "]";

    return result;
  }

  std::shared_ptr<Model> model(const std::string &type) {
    if (m_models.find(type) == m_models.end()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s model \"%s\" is not available.\n"
                                    "Available models:    %s\n",
                                    m_option.c_str(), type.c_str(),
                                    key_list(m_models).c_str());
    }

    return m_models[type]->create(m_grid);
  }

  template<class T>
  std::shared_ptr<Model> modifier(const std::string &type, std::shared_ptr<T> input) {
    if (m_modifiers.find(type) == m_modifiers.end()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s modifier \"%s\" is not available.\n"
                                    "Available modifiers:    %s\n",
                                    m_option.c_str(), type.c_str(),
                                    key_list(m_modifiers).c_str());
    }

    return m_modifiers[type]->create(m_grid, input);
  }

  // virtual base class that allows storing different model creators
  // in the same dictionary
  class ModelCreator {
  public:
    virtual std::shared_ptr<Model> create(IceGrid::ConstPtr g) = 0;
    virtual ~ModelCreator() {}
  };

  // Creator for a specific model class M.
  template <class M>
  class SpecificModelCreator : public ModelCreator {
  public:
    std::shared_ptr<Model> create(IceGrid::ConstPtr g) {
      return std::shared_ptr<Model>(new M(g));
    }
  };

  // virtual base class that allows storing different modifier
  // creators in the same dictionary
  class ModifierCreator {
  public:
    virtual std::shared_ptr<Model> create(IceGrid::ConstPtr g, std::shared_ptr<Model> input) = 0;
    virtual ~ModifierCreator() {}
  };

  // Creator for a specific modifier class M.
  template <class M>
  class SpecificModifierCreator : public ModifierCreator {
  public:
    std::shared_ptr<Model> create(IceGrid::ConstPtr g, std::shared_ptr<Model> input) {
      return std::shared_ptr<Model>(new M(g, input));
    }
  };

  std::string m_default_type, m_option;
  std::map<std::string, std::shared_ptr<ModelCreator> > m_models;
  std::map<std::string, std::shared_ptr<ModifierCreator> > m_modifiers;
  IceGrid::ConstPtr m_grid;
};

} // end of namespace pism

#endif /* _PCFACTORY_H_ */
