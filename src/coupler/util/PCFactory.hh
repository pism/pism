// Copyright (C) 2011, 2013, 2014 PISM Authors
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

#include "pism_const.hh"
#include "pism_options.hh"

#include "IceGrid.hh"
#include <map>

#include "error_handling.hh"

namespace pism {

class Config;

template <class Model, class Modifier>
class PCFactory {
public:

  // virtual base class that allows storing different model creators
  // in the same dictionary
  class ModelCreator {
  public:
    virtual Model* create(IceGrid &g) = 0;
    virtual ~ModelCreator() {}
  };

  // Creator for a specific model class M.
  template <class M>
  class SpecificModelCreator : public ModelCreator {
  public:
    M* create(IceGrid &g) {
      return new M(g);
    }
  };

  // virtual base class that allows storing different modifier
  // creators in the same dictionary
  class ModifierCreator {
  public:
    virtual Modifier* create(IceGrid &g, Model* input) = 0;
    virtual ~ModifierCreator() {}
  };

  // Creator for a specific modifier class M.
  template <class M>
  class SpecificModifierCreator : public ModifierCreator {
  public:
    M* create(IceGrid &g, Model* input) {
      return new M(g, input);
    }
  };

#ifdef PISM_USE_TR1
  typedef std::tr1::shared_ptr<ModelCreator> ModelCreatorPtr;
  typedef std::tr1::shared_ptr<ModifierCreator> ModifierCreatorPtr;
#else
  typedef std::shared_ptr<ModelCreator> ModelCreatorPtr;
  typedef std::shared_ptr<ModifierCreator> ModifierCreatorPtr;
#endif

  PCFactory<Model,Modifier>(IceGrid &g)
  : m_grid(g) {}
  ~PCFactory<Model,Modifier>() {}

  //! Sets the default type name.
  void set_default(std::string name) {
    if (m_models.find(name) == m_models.end()) {
      throw RuntimeError::formatted("type %s is not registered", name.c_str());
    } else {
      m_default_type = name;
    }
  }

  //! Creates a boundary model. Processes command-line options.
  Model* create() {
    std::vector<std::string> choices;
    std::string model_list, modifier_list, descr;
    bool flag = false;
    Model* result = NULL;

    // build a list of available models:
    typename std::map<std::string, ModelCreatorPtr >::iterator k;
    k = m_models.begin();
    model_list = "[" + (k++)->first;
    for (; k != m_models.end(); k++) {
      model_list += ", " + k->first;
    }
    model_list += "]";

    // build a list of available modifiers:
    typename std::map<std::string, ModifierCreatorPtr >::iterator p;
    p = m_modifiers.begin();
    modifier_list = "[" + (p++)->first;
    for (; p != m_modifiers.end(); p++) {
      modifier_list += ", " + p->first;
    }
    modifier_list += "]";

    descr =  "Sets up the PISM " + m_option + " model. Available models: " + model_list +
      " Available modifiers: " + modifier_list;

    // Get the command-line option:
    OptionsStringArray("-" + m_option, descr, m_default_type, choices, flag);

    if (choices.empty()) {
      if (flag == true) {
        throw RuntimeError::formatted("option -%s requires an argument.\n", m_option.c_str());
      }
      choices.push_back(m_default_type);
    }

    // the first element has to be an *actual* model (not a modifier), so we
    // create it:
    std::vector<std::string>::iterator j = choices.begin();

    if (m_models.find(*j) == m_models.end()) {
      throw RuntimeError::formatted("%s model \"%s\" is not available.\n"
                                    "Available models:    %s\n"
                                    "Available modifiers: %s",
                                    m_option.c_str(), j->c_str(),
                                    model_list.c_str(), modifier_list.c_str());
    }

    result = m_models[*j]->create(m_grid);

    ++j;

    // process remaining arguments:
    while (j != choices.end()) {
      if (m_modifiers.find(*j) == m_modifiers.end()) {
        throw RuntimeError::formatted("%s modifier \"%s\" is not available.\n"
                                      "Available modifiers: %s",
                                      m_option.c_str(), j->c_str(), modifier_list.c_str());
      }

      result =  m_modifiers[*j]->create(m_grid, result);

      ++j;
    }

    return result;
  }

  //! Adds a boundary model to the dictionary.
  template <class M>
  void add_model(std::string name) {
    m_models[name] = ModelCreatorPtr(new SpecificModelCreator<M>);
  }

  template <class M>
  void add_modifier(std::string name) {
    m_modifiers[name] = ModifierCreatorPtr(new SpecificModifierCreator<M>);
  }

  //! Removes a boundary model from the dictionary.
  void remove_model(std::string name) {
    m_models.erase(name);
  }

  void remove_modifier(std::string name) {
    m_modifiers.erase(name);
  }

  //! Clears the dictionary.
  void clear_models() {
    m_models.clear();
  }

  void clear_modifiers() {
    m_modifiers.clear();
  }
protected:
  virtual void add_standard_types() {}
  std::string m_default_type, m_option;
  std::map<std::string, ModelCreatorPtr> m_models;
  std::map<std::string, ModifierCreatorPtr> m_modifiers;
  IceGrid &m_grid;
};

} // end of namespace pism

#endif /* _PCFACTORY_H_ */
