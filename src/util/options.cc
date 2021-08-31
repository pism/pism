/* Copyright (C) 2014, 2015, 2016, 2017, 2018, 2020, 2021 PISM Authors
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

#include <cstring> // memset

#include <petscsys.h>

#include "error_handling.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Units.hh"
#include "pism/util/pism_utilities.hh"
#include "pism_options.hh"

namespace pism {
namespace options {

String::String(const std::string& option,
               const std::string& description) {
  int errcode = process(option, description, "", DONT_ALLOW_EMPTY);
  if (errcode != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "failed to process option %s", option.c_str());
  }
}

String::String(const std::string& option,
               const std::string& description,
               const std::string& default_value,
               ArgumentFlag argument_flag) {
  int errcode = process(option, description, default_value, argument_flag);
  if (errcode != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "failed to process option %s", option.c_str());
  }
}

static const int TEMPORARY_STRING_LENGTH = 32768;

int String::process(const std::string& option,
                    const std::string& description,
                    const std::string& default_value,
                    ArgumentFlag argument_flag) {

  char string[TEMPORARY_STRING_LENGTH];
  memset(string, 0, TEMPORARY_STRING_LENGTH);

  PetscBool flag = PETSC_FALSE;

  PetscErrorCode ierr;
  ierr = PetscOptionsGetString(NULL, // default option database
                               NULL, // no prefix
                               option.c_str(),
                               string,
                               TEMPORARY_STRING_LENGTH,
                               &flag);
  PISM_CHK(ierr, "PetscOptionsGetString");

  std::string result = string;

  if (flag == PETSC_TRUE) {
    if (result.empty()) {
      if (argument_flag == ALLOW_EMPTY) {
        this->set("", true);
      } else {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "command line option '%s'\n"
                                      "(%s)\n"
                                      "requires an argument.",
                                      option.c_str(), description.c_str());
      }
    } else {
      this->set(result, true);
    }
  } else {
    this->set(default_value, false);
  }

  return 0;
}

Keyword::Keyword(const std::string& option,
                 const std::string& description,
                 const std::string& choices,
                 const std::string& default_value) {

  if (choices.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "empty choices argument");
  }

  std::string list = "[" + choices + "]";
  std::string long_description = description + " Choose one of " + list;

  String input(option, long_description, default_value, DONT_ALLOW_EMPTY);

  // use the default value if the option was not set
  if (not input.is_set()) {
    this->set(input, input.is_set());
    return;
  }

  std::string word = input;
  // find ":" and discard everything that goes after
  size_t n = word.find(":");
  if (n != std::string::npos) {
    word.resize(n);
  }

  // transform a comma-separated list of choices into a set of
  // choices:
  auto choices_set = set_split(choices, ',');

  // use the choice if it is valid and stop if it is not
  if (choices_set.find(word) != choices_set.end()) {
    this->set(word, true);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid %s argument: '%s'. Please choose one of %s.\n",
                                  option.c_str(), word.c_str(), list.c_str());
  }
}

Integer::Integer(const std::string& option,
                 const std::string& description,
                 int default_value) {

  String input(option, description,
               pism::printf("%d", default_value),
               DONT_ALLOW_EMPTY);

  if (input.is_set()) {
    long int result = 0;
    try {
      result = parse_integer(input);
    } catch (RuntimeError &e) {
      e.add_context("processing command-line option '%s %s'",
                    option.c_str(), input->c_str());
      throw;
    }
    this->set(static_cast<int>(result), true);
  } else {
    this->set(static_cast<int>(default_value), false);
  }
}


IntegerList::IntegerList(const std::string& option,
                         const std::string& description,
                         const std::vector<int> &defaults) {
  std::vector<double> default_value;

  for (auto v : defaults) {
    default_value.push_back(v);
  }

  RealList input(option, description, default_value);
  std::vector<int> result;

  for (auto v : input.value()) {
    if (fabs(v - floor(v)) > 1e-6) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "Can't process '%s': (%f is not an integer).",
                                    option.c_str(), v);
    }
    result.push_back(static_cast<int>(v));
  }

  this->set(result, input.is_set());
}

const int& IntegerList::operator[](size_t index) const {
  return m_value[index];
}

Real::Real(std::shared_ptr<units::System> system,
           const std::string& option,
           const std::string& description,
           const std::string& units,
           double default_value) {

  std::string buffer = pism::printf("%f", default_value);

  String input(option, description, buffer, DONT_ALLOW_EMPTY);

  if (input.is_set()) {
    char *endptr = NULL;
    double result = strtod(input->c_str(), &endptr);
    if (*endptr != '\0') {
      // assume that "input" contains units and try converting to "units":
      try {
        result = units::convert(system, 1.0, input.value(), units);
      } catch (RuntimeError &e) {
        e.add_context("trying to convert '%s' to '%s'",
                      input->c_str(), units.c_str());
        e.add_context("processing the command-line option %s",
                      option.c_str());
        throw;
      }
    }
    this->set(result, true);
  } else {
    this->set(default_value, false);
  }
}


RealList::RealList(const std::string& option,
                   const std::string& description,
                   const std::vector<double> &default_value) {
  String input(option, description, "", DONT_ALLOW_EMPTY);
  std::vector<double> result = default_value;

  if (input.is_set()) {
    result.clear();
    for (auto p : split(input, ',')) {
      result.push_back(parse_number(p));
    }
  }
  this->set(result, input.is_set());
}

const double& RealList::operator[](size_t index) const {
  return m_value[index];
}

bool Bool(const std::string& option,
          const std::string& description) {
  return String(option, description, "", ALLOW_EMPTY).is_set();
}

//! Stop if an option `old_name` is set, printing a message that `new_name` should be used instead.
void deprecated(const std::string &old_name, const std::string &new_name) {

  String option(old_name, "no description", "default",
                options::ALLOW_EMPTY);

  if (option.is_set()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "command-line option '%s' is deprecated."
                                  " Please use '%s' instead.",
                                  old_name.c_str(), new_name.c_str());
  }
}

//! Print a warning telling the user that an option was ignored.
void ignored(const Logger &log, const std::string &name) {

  String option(name, "no description", "default");

  if (option.is_set()) {
    log.message(1, "PISM WARNING: ignoring command-line option '%s'.\n",
                name.c_str());
  }
}

//!Stop if an option `name` is set.
void forbidden(const std::string &name) {
  bool option_is_set = options::Bool(name, "no description");

  if (option_is_set) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "command-line option '%s' is not allowed.",
                                  name.c_str());
  }
}

} // end of namespace options
} // end of namespace pism
