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
                         const std::string &IOmode,
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
  int mode    = define_mode(IOmode);
  m_comp_comm = pioInit(glob, n_writers, mode,
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

int Initializer::define_mode(const std::string &mode) {
  std::map<std::string, int> modes =
    {
     {"PIO_NONE", PIO_NONE},
     {"PIO_MPI", PIO_MPI},
     {"PIO_WRITER", PIO_WRITER},
     {"PIO_ASYNCH", PIO_ASYNCH},
     {"PIO_FPGUARD", PIO_FPGUARD},
     {"PIO_MPI_FW_ORDERED", PIO_MPI_FW_ORDERED},
     {"PIO_MPI_FW_AT_ALL", PIO_MPI_FW_AT_ALL},
     {"PIO_MPI_FW_AT_REBLOCK", PIO_MPI_FW_AT_REBLOCK}
    };

  if (modes.find(mode) == modes.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid CDI-PIO I/O mode: %s", mode.c_str());
  }

  return modes[mode];
}

} // namespace cdipio
} // namespace pism
