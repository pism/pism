// Copyright (C) 2011, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2024, 2025 PISM Authors
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

#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Config.hh"
#include "pism/util/Context.hh"

namespace pism {

template <class Model>
class PCFactory {
public:

  PCFactory(std::shared_ptr<const Grid> g, const std::string &parameter)
    : m_parameter(parameter), m_grid(g)  {}
  ~PCFactory() {}

  //! Creates a boundary model. Processes command-line options.
  virtual std::shared_ptr<Model> create() {

    auto choices = m_grid->ctx()->config()->get_string(m_parameter);

    return create(choices);
  }

  void validate(const std::string &list) const {

    auto choices = split(list, ',');

    auto doc = m_grid->ctx()->config()->doc(m_parameter);
    auto opt = m_grid->ctx()->config()->option(m_parameter);

    if (choices.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "Parameter %s (%s) cannot be empty (got '%s')",
                                    m_parameter.c_str(), doc.c_str(), list.c_str());
    }

    auto model1 = m_models.begin()->first;

    if (m_models.find(choices[0]) == m_models.end() and
        m_modifiers.find(choices[0]) != m_modifiers.end()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "The first item on the list %s\n"
                                    "(%s)\n"
                                    "has to be a 'model' (one of %s),\n"
                                    "while the rest have to be 'modifiers' (one of %s).\n"
                                    "Got %s='%s'.\n"
                                    "To use %s you also have to select a model, e.g. using the command-line option\n"
                                    "'-%s %s,%s'.",
                                    m_parameter.c_str(), doc.c_str(), key_list(m_models).c_str(),
                                    key_list(m_modifiers).c_str(),
                                    m_parameter.c_str(), list.c_str(), list.c_str(),
                                    opt.c_str(), model1.c_str(), list.c_str());
    }
  }

  //! Creates a boundary model.
  virtual std::shared_ptr<Model> create(const std::string &type) {
    validate(type);

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

  template <typename T>
  std::string key_list(std::map<std::string, T> list) const {
    std::vector<std::string> keys;

    for (const auto &i : list) {
      keys.push_back(i.first);
    }

    return "[" + join(keys, ", ") + "]";
  }

  std::shared_ptr<Model> model(const std::string &type) {
    if (m_models.find(type) == m_models.end()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "cannot allocate %s \"%s\".\n"
                                    "Available models:    %s\n",
                                    m_parameter.c_str(), type.c_str(), key_list(m_models).c_str());
    }

    return m_models[type]->create(m_grid);
  }

  template <class T>
  std::shared_ptr<Model> modifier(const std::string &type, std::shared_ptr<T> input) {
    if (m_modifiers.find(type) == m_modifiers.end()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "cannot allocate %s modifier \"%s\".\n"
                                    "Available modifiers:    %s\n",
                                    m_parameter.c_str(), type.c_str(),
                                    key_list(m_modifiers).c_str());
    }

    return m_modifiers[type]->create(m_grid, input);
  }

  // virtual base class that allows storing different model creators
  // in the same dictionary
  class ModelCreator {
  public:
    virtual std::shared_ptr<Model> create(std::shared_ptr<const Grid> g) = 0;
    virtual ~ModelCreator() = default;
  };

  // Creator for a specific model class M.
  template <class M>
  class SpecificModelCreator : public ModelCreator {
  public:
    std::shared_ptr<Model> create(std::shared_ptr<const Grid> g) {
      return std::shared_ptr<Model>(new M(g));
    }
  };

  // virtual base class that allows storing different modifier
  // creators in the same dictionary
  class ModifierCreator {
  public:
    virtual std::shared_ptr<Model> create(std::shared_ptr<const Grid> g,
                                          std::shared_ptr<Model> input) = 0;

    virtual ~ModifierCreator() = default;
  };

  // Creator for a specific modifier class M.
  template <class M>
  class SpecificModifierCreator : public ModifierCreator {
  public:
    std::shared_ptr<Model> create(std::shared_ptr<const Grid> g, std::shared_ptr<Model> input) {
      return std::shared_ptr<Model>(new M(g, input));
    }
  };

  std::string m_parameter;
  std::map<std::string, std::shared_ptr<ModelCreator> > m_models;
  std::map<std::string, std::shared_ptr<ModifierCreator> > m_modifiers;
  std::shared_ptr<const Grid> m_grid;
};

} // end of namespace pism

#endif /* _PCFACTORY_H_ */
