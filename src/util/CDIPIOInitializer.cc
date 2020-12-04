#include <mpi.h>
#if (Pism_USE_CDIPIO==1)
extern "C"{
#include "cdipio.h"
#include "cdi.h"
}
#endif



namespace pism {
namespace cdipio {

Initializer::Initializer(int n_writers, int IOmode, MPI_Comm glob) {
#if (Pism_USE_CDIPIO==1)
	// Initialize YAXT library
	xt_initialize(glob);
	// Initialize CDI-PIO library
	float partInflate = 1.0;
	m_local_comm = pioInit(glob, n_writers, IOmode, &m_pioNamespace, partInflate, cdiPioNoPostCommSetup);
#endif
}

Initializer::~Initializer() {
#if (Pism_USE_CDIPIO==1)
	// Finalize CDI-PIO library
	pioFinalize();
	// Finalize YAXT library
	xt_finalize();
#endif
}

MPI_Comm Initializer::get_comp_comm() {
#if (Pism_USE_CDIPIO==1)
	return m_local_comm;
#endif
}

int Initializer::get_pionamespace() {
#if (Pism_USE_CDIPIO==1)
	return m_pioNamespace;
#endif
}

void Initializer::activate_namespace() {
	namespaceSetActive(m_pioNamespace);
}

}
}