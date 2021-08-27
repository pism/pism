/* Copyright (C) 2021 PISM Authors
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
#ifndef CDIPIOWINDOW_H
#define CDIPIOWINDOW_H

#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace cdipio {

class Window {
public:
  Window(IO_Backend b);
  void update(bool sth);
  void expose();
private:
  bool m_sthwritten;
  IO_Backend m_backend;
};

} // end of namespace cdipio
} // end of namespace pism

#endif /* CDIPIOWINDOW */
