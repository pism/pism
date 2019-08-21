/* Copyright (C) 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include <sstream>
#include <petscsys.h>

#include "pism_options.hh"
#include "error_handling.hh"
#include "pism/util/Logger.hh"
#include "pism/util/pism_utilities.hh"

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

  char tmp[TEMPORARY_STRING_LENGTH];
  PetscBool flag = PETSC_FALSE;

  memset(tmp, 0, TEMPORARY_STRING_LENGTH);

  PetscErrorCode ierr;
  ierr = PetscOptionsBegin(MPI_COMM_SELF, "", "", "");
  PISM_CHK(ierr, "PetscOptionsBegin");

  ierr = PetscOptionsString(option.c_str(),
                            description.c_str(),
                            "", // manual page
                            default_value.c_str(), // default value
                            tmp,                   // output
                            TEMPORARY_STRING_LENGTH, // max. length of the output
                            &flag);                  // PETSC_TRUE if found, else PETSC_FALSE
  PISM_CHK(ierr, "PetscOptionsString");

  ierr = PetscOptionsEnd();
  PISM_CHK(ierr, "PetscOptionsEnd");

  std::string result = tmp;

  if (flag == PETSC_TRUE) {
    if (result.empty()) {
      if (argument_flag == ALLOW_EMPTY) {
        this->set("", true);
      } else {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "command line option '%s' requires an argument.",
                                      option.c_str());
      }
    } else {
      this->set(result, true);
    }
  } else {
    this->set(default_value, false);
  }

  return 0;
}

StringList::StringList(const std::string& option,
                       const std::string& description,
                       const std::string& default_value) {
  String input(option, description, default_value, DONT_ALLOW_EMPTY);

  this->set(split(input, ','), input.is_set());
}

const std::string& StringList::operator[](size_t index) const {
  return m_value[index];
}

std::string StringList::to_string() {
  return join(m_value, ",");
}

StringSet::StringSet(const std::string& option,
                     const std::string& description,
                     const std::string& default_value) {
  StringList input(option, description, default_value);
  std::set<std::string> result;

  for (auto s : input.value()) {
    result.insert(s);
  }

  this->set(result, input.is_set());
}

std::string StringSet::to_string() {
  return set_join(m_value, ",");
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
  Real input(option, description, default_value);
  double result = input;
  if (fabs(result - floor(result)) > 1e-6) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't process '%s': (%f is not an integer).",
                                  option.c_str(), result);
  }
  this->set(static_cast<int>(result), input.is_set());
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

Real::Real(const std::string& option,
           const std::string& description,
           double default_value) {

  std::stringstream buffer;
  // NB! This may round default_value.
  buffer << default_value;

  String input(option, description, buffer.str(), DONT_ALLOW_EMPTY);

  if (input.is_set()) {
    char *endptr = NULL;
    double result = strtod(input->c_str(), &endptr);
    if (*endptr != '\0') {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't parse '%s %s': (%s is not a number).",
                                    option.c_str(), input->c_str(), input->c_str());
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
    std::istringstream arg(input);
    std::string tmp;

    result.clear();
    while(getline(arg, tmp, ',')) {
      double d;
      char *endptr;

      d = strtod(tmp.c_str(), &endptr);
      if (*endptr != '\0') {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't parse %s (%s is not a number).",
                                      tmp.c_str(), tmp.c_str());
      } else {
        result.push_back(d);
      }
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
