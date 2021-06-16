#include "CDIPIOInitializer.hh"

#include <mpi.h>

extern "C" {
#include <cdi.h>
#include <cdipio.h>
#include <yaxt.h>
}

namespace pism {
namespace cdipio {

Initializer::Initializer(int n_writers, int IOmode, MPI_Comm glob) {
  // check if CDIPIO is used
  if (n_writers > 0) {
    m_async = true;
  } else {
    m_async = false;
  }
  
  // initialize CDIPIO library
  if (m_async) {
    xt_initialize(glob);
    float partInflate = 1.0;
    m_comp_comm      = pioInit(glob, n_writers, IOmode,
                               &m_pioNamespace, partInflate, cdiPioNoPostCommSetup);
  } else {
    m_comp_comm    = MPI_COMM_WORLD;
    m_pioNamespace = -1;
  }
}

Initializer::~Initializer() {
  // finalize CDIPIO and YAXT libraries
  if (m_async) {
    pioFinalize();
    xt_finalize();
  }
}

MPI_Comm Initializer::comp_comm() {
  return m_comp_comm;
}

int Initializer::pio_namespace() {
  return m_pioNamespace;
}

void Initializer::activate_namespace() {
  if (m_async) {
    namespaceSetActive(m_pioNamespace);
  }
}

} // namespace cdipio
} // namespace pism
