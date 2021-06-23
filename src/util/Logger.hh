/* Copyright (C) 2015, 2016, 2021 PISM Authors
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
#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <string>
#include <memory>

#include <mpi.h>                // MPI_Comm

namespace pism {

typedef enum {WARNING=1, DEBUG2=2, DEBUG3=3, TRACE=4} LoggerLevel;

//! A basic logging class.
/**
 * The default implementation (message_impl()) just prints to `stdout` on rank 0 of the
 * communicator.
 *
 * This class was created to make it possible to silence PISM's output when it is used as a library
 * and make it possible to separate outputs from different PISM (IceModel, etc) instances running
 * side by side.
 */
class Logger {
public:
  Logger(MPI_Comm com, int threshold);
  virtual ~Logger();

  typedef std::shared_ptr<Logger> Ptr;
  typedef std::shared_ptr<const Logger> ConstPtr;

  //! Print a message to the log.
  /** Does nothing if `threshold` is greater than the value provided to the constructor or set using
   *  set_threshold().
   */
  void message(int threshold, const char format[], ...) const __attribute__((format(printf, 3, 4)));
  void message(int threshold, const std::string &text) const;

  //! Print an error message to the log.
  /** Always prints the message (regardless of the threshold). The base class implementation prints
   * to stderr.
   */
  void error(const char format[], ...) const __attribute__((format(printf, 2, 3)));

  //! Set verbosity threshold.
  void set_threshold(int level);

  //! Get verbosity threshold.
  int get_threshold() const;

  //! Silence the logger.
  /**
   * This makes it possible to temporarily silence the logger and then re-enable it with the same
   * threshold as before, but without explicitly storing the current threshold.
   */
  void disable() const;
  //! (Re-)enable the logger.
  void enable() const;
protected:
  //! Do the hard work. Override this in a derived class to customize.
  virtual void message_impl(const char buffer[]) const;
  virtual void error_impl(const char buffer[]) const;
  private:
  struct Impl;
  Impl *m_impl;
  Logger(const Logger&);
  Logger & operator=(const Logger &);
};

//! A logger that accumulates messages and reports them as a string.
class StringLogger : public Logger {
public:
  StringLogger(MPI_Comm com, int threshold);
  virtual ~StringLogger();

  void reset();

  std::string get() const;
protected:
  virtual void message_impl(const char buffer[]) const;
  virtual void error_impl(const char buffer[]) const;
private:
  struct Impl;
  Impl *m_impl;
};

Logger::Ptr logger_from_options(MPI_Comm com);

} // end of namespace pism

#endif /* _LOGGER_H_ */
