#include "CDIPIOInitializer.hh"

#include "pism/util/IceGrid.hh"

#include "pism/util/error_handling.hh"

#include <mpi.h>
#if (Pism_USE_CDIPIO==1)
extern "C"{
#include "cdipio.h"
#include "cdi.h"
}
#endif



namespace pism {
namespace cdipio {

Initializer::Initializer(int n_writers, int IOmode, MPI_Comm glob, bool async) {
#if (Pism_USE_CDIPIO==1)
	m_async = async;
	if (m_async) {
		xt_initialize(glob);
		float partInflate = 1.0;
		m_local_comm = pioInit(glob, n_writers, IOmode, &m_pioNamespace, partInflate, cdiPioNoPostCommSetup);
	} else {
		m_local_comm = MPI_COMM_WORLD;
		m_pioNamespace = -1;
	}
#endif
}

Initializer::~Initializer() {
#if (Pism_USE_CDIPIO==1)
	if (m_async) {
		pioFinalize();
		xt_finalize();
	}
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
	if (m_async) namespaceSetActive(m_pioNamespace);
}

}
}
