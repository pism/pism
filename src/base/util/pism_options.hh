// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#ifndef _PISM_OPTIONS_H_
#define _PISM_OPTIONS_H_

#include "pism_const.hh"

namespace pism {

class Config;

namespace options {

template <typename T>
class Option {
public:
  Option() {
    m_is_set = false;
  }
  Option(T value, bool is_set)
    : m_value(value) {
    m_is_set = is_set;
  }
  operator T() {
    return m_value;
  }
  bool is_set() {
    return m_is_set;
  }
  T value() {
    return m_value;
  }
private:
  T m_value;
  bool m_is_set;
};

typedef Option<std::string> String;
typedef Option<std::vector<std::string> > StringList;
typedef Option<std::set<std::string> > StringSet;

typedef Option<int> Integer;
typedef Option<std::vector<int> > IntegerList;

typedef Option<double> Real;
typedef Option<std::vector<double> > RealList;

String string(const std::string& option,
              const std::string& description,
              const std::string& default_value,
              bool allow_empty_arg = false);

StringList string_list(const std::string& option,
                       const std::string& description,
                       const std::string& default_value);

StringSet string_set(const std::string& option,
                     const std::string& description,
                     const std::string& default_value);

String keyword(const std::string& option,
               const std::string& description,
               const std::string& choices,
               const std::string& default_value);

Integer integer(const std::string& option,
                const std::string& description,
                int default_value);

IntegerList integer_list(const std::string& option,
                         const std::string& description);

Real real(const std::string& option,
          const std::string& description,
          double default_value);

RealList real_list(const std::string& option,
                   const std::string& description);

} // end of namespace options

PetscErrorCode verbosityLevelFromOptions();

// handy functions for processing options:
PetscErrorCode OptionsList(std::string opt, std::string text, std::set<std::string> choices,
                           std::string default_value, std::string &result, bool &flag);

PetscErrorCode OptionsString(std::string option, std::string text,
                             std::string &result, bool &flag, bool allow_empty_arg = false);
PetscErrorCode OptionsStringArray(std::string opt, std::string text, std::string default_value,
                                      std::vector<std::string>& result, bool &flag);
PetscErrorCode OptionsStringSet(std::string opt, std::string text, std::string default_value,
                                std::set<std::string>& result, bool &flag);

PetscErrorCode OptionsInt(std::string option, std::string text,
                              int &result, bool &is_set);
PetscErrorCode OptionsIntArray(std::string option, std::string text,
                                   std::vector<int> &result, bool &is_set);

PetscErrorCode OptionsReal(std::string option, std::string text,
                               double &result, bool &is_set);
PetscErrorCode OptionsRealArray(std::string option, std::string text,
                                    std::vector<double> &result, bool &is_set);

bool OptionsIsSet(std::string option);
bool OptionsIsSet(std::string option, std::string descr);

PetscErrorCode OptionsHasArgument(std::string option, bool &result);

PetscErrorCode ignore_option(MPI_Comm com, std::string name);
PetscErrorCode check_old_option_and_stop(std::string old_name, std::string new_name);
PetscErrorCode stop_if_set(std::string name);

// usage message and required options; drivers use these
PetscErrorCode stop_on_version_option();

PetscErrorCode show_usage_and_quit(MPI_Comm com, std::string execname, std::string usage);

PetscErrorCode show_usage_check_req_opts(MPI_Comm com, std::string execname,
                                         std::vector<std::string> required_options,
                                         std::string usage);

// config file initialization:
PetscErrorCode init_config(MPI_Comm com,
                           Config &config, Config &overrides,
                           bool process_options = false);

PetscErrorCode set_config_from_options(Config &config);


} // end of namespace pism

#endif /* _PISM_OPTIONS_H_ */
