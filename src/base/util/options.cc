/* Copyright (C) 2014 PISM Authors
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

#include "pism_options.hh"
#include "error_handling.hh"

namespace pism {
namespace options {

String string(const std::string& option,
              const std::string& description,
              const std::string& default_value,
              bool allow_empty_arg) {

  char tmp[TEMPORARY_STRING_LENGTH];
  PetscBool flag = PETSC_FALSE;

  PetscErrorCode ierr;
  ierr = PetscOptionsString(option.c_str(),
                            description.c_str(),
                            "", // manual page
                            default_value.c_str(), // default value
                            tmp,                   // output
                            TEMPORARY_STRING_LENGTH, // max. length of the output
                            &flag);                  // PETSC_TRUE if found, else PETSC_FALSE
  PISM_PETSC_CHK(ierr, "PetscOptionsString");

  std::string result = tmp;

  if (flag == PETSC_TRUE) {
    if (result.empty()) {
      if (allow_empty_arg) {
        return String("", true);
      } else {
        throw RuntimeError::formatted("command line option '%s' requires an argument.",
                                      option.c_str());
      }
    } else {
      return String(result, true);
    }
  } else {
    return String(default_value, false);
  }
}


StringList string_list(const std::string& option,
                       const std::string& description,
                       const std::string& default_value) {
  String input = string(option, description, default_value, false);
  std::vector<std::string> result;

  std::string token;
  std::istringstream arg(input);
  while (getline(arg, token, ',')) {
    result.push_back(token);
  }

  return StringList(result, input.is_set());
}


StringSet string_set(const std::string& option,
                     const std::string& description,
                     const std::string& default_value) {
  StringList input = string_list(option, description, default_value);
  std::set<std::string> result;

  std::vector<std::string> array = input;
  for (unsigned int k = 0; k < array.size(); ++k) {
    result.insert(array[k]);
  }

  return StringSet(result, input.is_set());
}

String keyword(const std::string& option,
               const std::string& description,
               const std::string& choices,
               const std::string& default_value) {

  if (choices.empty()) {
    throw RuntimeError::formatted("empty choices argument");
  }

  std::string list = "[" + choices + "]";
  std::string long_description = description + " Choose one of " + list;

  String input = string(option, long_description, default_value, false);
  
  // return the default value if the option was not set
  if (not input.is_set()) {
    return input;
  }

  std::string word = input;
  // find ":" and discard everything that goes after
  size_t n = word.find(":");
  if (n != std::string::npos) {
    word.resize(n);
  }

  // transform a comma-separated list of choices into a set of
  // choices:
  std::set<std::string> choices_set;
  {  
    std::string token;
    std::istringstream arg(choices);
    while (getline(arg, token, ',')) {
      choices_set.insert(token);
    }
  }

  // return the choice if it is valid and stop if it is not
  if (choices_set.find(word) != choices_set.end()) {
    return String(word, true);
  } else {
    throw RuntimeError::formatted("invalid %s argument: '%s'. Please choose one of %s.\n",
                                  option.c_str(), word.c_str(), list.c_str());
  }
}




Integer integer(const std::string& option,
                const std::string& description,
                int default_value) {
  Real input = real(option, description, default_value);
  double result = input;
  if (fabs(result - floor(result)) > 1e-6) {
    throw RuntimeError::formatted("Can't process '%s': (%f is not an integer).",
                                  option.c_str(), result);
  }
  return Integer((int)result, input.is_set());
}


IntegerList integer_list(const std::string& option,
                         const std::string& description) {
  RealList input = real_list(option, description);
  std::vector<int> result;
  std::vector<double> array = input;

  for (unsigned int k = 0; k < array.size(); ++k) {
    if (fabs(array[k] - floor(array[k])) > 1e-6) {
      throw RuntimeError::formatted("Can't process '%s': (%f is not an integer).",
                                    option.c_str(), array[k]);
    }
    result.push_back((int)array[k]);
  }
  
  return IntegerList(result, input.is_set());
}


Real real(const std::string& option,
          const std::string& description,
          double default_value) {
  std::stringstream buffer;
  buffer << default_value;
  String input = string(option, description, buffer.str(), false);

  std::string str = input;
  char *endptr = NULL;
  double result = strtod(str.c_str(), &endptr);
  if (*endptr != '\0') {
    throw RuntimeError::formatted("Can't parse '%s %s': (%s is not a number).",
                                  option.c_str(), str.c_str(), str.c_str());
  }
  return Real(result, input.is_set());
}


RealList real_list(const std::string& option,
                   const std::string& description) {
  String input = string(option, description, "", false);
  std::vector<double> result;

  if (input.is_set()) {
    std::istringstream arg(input);
    std::string tmp;

    result.clear();
    while(getline(arg, tmp, ',')) {
      double d;
      char *endptr;

      d = strtod(tmp.c_str(), &endptr);
      if (*endptr != '\0') {
        throw RuntimeError::formatted("Can't parse %s (%s is not a number).",
                                      tmp.c_str(), tmp.c_str());
      } else {
        result.push_back(d);
      }
    }
  }
  return RealList(result, input.is_set());
}

} // end of namespace options
} // end of namespace pism
