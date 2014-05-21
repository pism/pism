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

#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <string>

namespace pism {

/** Database of configuration flags and parameters for PISM.
 *
 * Stores scalar (type 'double'), boolean, and string parameters.
 *
 * Requirements:
 *
 * 1. Self-documenting. We should be able to generate documentation
 * using the database that describes the logic. Documentation should
 * allow including references to literature.
 *
 * 2. Support unit conversion for physical constants and parameters.
 *
 * 3. Support parameter studies by allowing a user to provide a subset
 * of the database with parameters relevant in a particular case.
 * ("Overrides")
 *
 * 4. Establish an one-to-one mapping from configuration parameters to
 * command-line options.
 *
 * 5. Allow short versions of command-line options.
 *
 * 6. The user should be able to disable processing of command-line
 * options. (Makes it easier to use as a library.)
 *
 * 7. Support command-line option prefixes to allow running more than
 * one instance.
 *
 * 8. Keep track of parameters that were set using command-line
 * options or an override file.
 *
 * 9. Provide a purely serial interface, avoid I/O, including standard
 * output and standard error.
 *
 */
class Config {
public:
  virtual ~Config();

  /** Set a scalar parameter `name` to the value `value`.
   *
   * Converts parameter from `units` to appropriate units used
   * internally.
   */
  void set_double(const std::string &name, double value,
                  const std::string &units);

  /** Set a Boolean flag `name` to `value`.
   */
  void set_boolean(const std::string &name, bool value);

  /** Set a string parameter `name` to `value`.
   */
  void set_string(const std::string &name, const std::string &value);

  /** Get a scalar parameter `name`, converting it to `units`.
   */
  double get_double(const std::string &name,
                    const std::string &units) const;

  /** Get a string parameter `name`.
   */
  std::string get_string(const std::string &name) const;

  /** Get a Boolean flag `name`.
   */
  bool get_boolean(const std::string &name) const;
private:
  virtual void set_double_impl(const std::string &name, double value) = 0;
  virtual void set_boolean_impl(const std::string &name, bool value) = 0;
  virtual void set_string_impl(const std::string &name, const std::string &value) = 0;

  virtual double get_double_impl(const std::string &name) const = 0;
  virtual std::string get_string_impl(const std::string &name) const = 0;
  virtual bool get_boolean_impl(const std::string &name) const = 0;
};

} // end of namespace pism

#endif /* _CONFIG_H_ */
