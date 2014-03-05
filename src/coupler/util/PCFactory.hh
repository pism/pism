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

class PISMConfig;

template <class Model, class Modifier>
class PCFactory {
public:
  PCFactory<Model,Modifier>(IceGrid& g, const PISMConfig& conf)
  : grid(g), config(conf) {}
  virtual ~PCFactory<Model,Modifier>() {}

  //! Sets the default type name.
  virtual PetscErrorCode set_default(std::string name) {
    void (*func) (IceGrid&, const PISMConfig&, Model*&);

    func = models[name];
    if (!func) {
      SETERRQ1(grid.com, 1,"ERROR: type %s is not registered", name.c_str());
    } else {
      default_type = name;
    }
    return 0;
  }

  //! Creates a boundary model. Processes command-line options.
  virtual PetscErrorCode create(Model* &result) {
    void (*F) (IceGrid&, const PISMConfig&, Model*&);
    PetscErrorCode ierr;
    std::vector<std::string> choices;
    std::string model_list, modifier_list, descr;
    bool flag = false;

    // build a list of available models:
    typename std::map<std::string,void(*)(IceGrid&, const PISMConfig&, Model*&)>::iterator k;
    k = models.begin();
    model_list = "[" + (k++)->first;
    for(; k != models.end(); k++) {
      model_list += ", " + k->first;
    }
    model_list += "]";

    // build a list of available modifiers:
    typename std::map<std::string,void(*)(IceGrid&, const PISMConfig&, Model*, Modifier*&)>::iterator p;
    p = modifiers.begin();
    modifier_list = "[" + (p++)->first;
    for(; p != modifiers.end(); p++) {
      modifier_list += ", " + p->first;
    }
    modifier_list += "]";

    descr =  "Sets up the PISM " + option + " model. Available models: " + model_list +
      " Available modifiers: " + modifier_list;

    // Get the command-line option:
    ierr = PISMOptionsStringArray("-" + option, descr, default_type, choices, flag); CHKERRQ(ierr);

    if (choices.empty()) {
      if (flag) {
        PetscPrintf(grid.com, "ERROR: option -%s requires an argument.\n", option.c_str());
        PISMEnd();
      }
      choices.push_back(default_type);
    }

    // the first element has to be an *actual* model (not a modifier), so we
    // create it:
    std::vector<std::string>::iterator j = choices.begin();

    F = models[*j];
    if (!F) {
      PetscPrintf(grid.com,
                  "ERROR: %s model \"%s\" is not available.\n"
                  "  Available models:    %s\n"
                  "  Available modifiers: %s\n",
                  option.c_str(), j->c_str(),
                  model_list.c_str(), modifier_list.c_str());
      PISMEnd();
    }

    (*F)(grid, config, result);

    ++j;

    // process remaining arguments:
    while (j != choices.end()) {
      void (*M) (IceGrid&, const PISMConfig&, Model*, Modifier*&);
      Modifier *mod;

      M = modifiers[*j];
      if (!M) {
        PetscPrintf(grid.com,
                    "ERROR: %s modifier \"%s\" is not available.\n"
                    "  Available modifiers: %s\n",
                    option.c_str(), j->c_str(), modifier_list.c_str());
        PISMEnd();
      }

      (*M)(grid, config, result, mod);

      result = mod;

      ++j;
    }

    return 0;
  }

  //! Adds a boundary model to the dictionary.
  virtual void add_model(std::string name, void(*func)(IceGrid&, const PISMConfig&, Model*&)) {
    models[name] = func;
  }

  virtual void add_modifier(std::string name, void(*func)(IceGrid&, const PISMConfig&, Model*, Modifier*&)) {
    modifiers[name] = func;
  }

  //! Removes a boundary model from the dictionary.
  virtual void remove_model(std::string name) {
    models.erase(name);
  }

  virtual void remove_modifier(std::string name) {
    modifiers.erase(name);
  }

  //! Clears the dictionary.
  virtual void clear_models() {
    models.clear();
  }

  virtual void clear_modifiers() {
    modifiers.clear();
  }
protected:
  virtual void add_standard_types() {}
  std::string default_type, option;
  std::map<std::string,void(*)(IceGrid&, const PISMConfig&, Model*&)> models;
  std::map<std::string,void(*)(IceGrid&, const PISMConfig&, Model*, Modifier*&)> modifiers;
  IceGrid& grid;
  const PISMConfig& config;
};

#endif /* _PCFACTORY_H_ */
