#include "CDIPIOInitializer.hh"

#include <mpi.h>

extern "C" {
#include <cdi.h>                // namespaceSetActive()
#include <cdipio.h>
#include <yaxt.h>               // xt_initialize(), xt_finalize()
}

#include <map>
#include "pism/util/error_handling.hh"

namespace pism {
namespace cdipio {

Initializer::Initializer(unsigned int n_writers,
                         MPI_Comm glob)
  : m_comp_comm(glob),
    m_pioNamespace(-1),
    m_initialized(false) {

  if (n_writers == 0) {
    return;
  }

  // initialize CDIPIO library
  xt_initialize(glob);
  float partInflate = 1.0;
  m_comp_comm = pioInit(glob, n_writers, PIO_MPI_FW_AT_ALL,
                        &m_pioNamespace, partInflate, cdiPioNoPostCommSetup);

  m_initialized = true;
}

Initializer::~Initializer() {
  // finalize CDIPIO and YAXT libraries
  if (m_initialized) {
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
  if (m_initialized) {
    namespaceSetActive(m_pioNamespace);
  }
}

} // namespace cdipio
} // namespace pism
