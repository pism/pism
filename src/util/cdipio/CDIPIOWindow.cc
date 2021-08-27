#include "CDIPIOWindow.hh"

#include <mpi.h>
extern "C" {
#include <cdipio.h>
}


namespace pism {
namespace cdipio {

Window::Window(IO_Backend b)
  : m_sthwritten(false),
    m_backend(b) {
}

void Window::update(bool sth) {
  m_sthwritten = sth;
}

void Window::expose() {
  if (m_backend == PISM_CDI) {
    if (m_sthwritten) {
      pioWriteTimestep();
      m_sthwritten = false;
    }
  } else {
    return;
  }
}

} // namespace cdipio
} // namespace pism
