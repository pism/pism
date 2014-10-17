#include "pism_python_signal.hh"
#include <signal.h>
#include <vector>
#include <errno.h>
#include "petsc.h"
#include "pism_const.hh"

#include "error_handling.hh"

using namespace pism;

SigInstaller::SigInstaller(int sig, void (*new_handler)(int) )
{
  m_old_handler = signal(sig, new_handler);
  m_sig = sig;
}

void SigInstaller::release()
{
  if (m_old_handler)
  {
    signal(m_sig, m_old_handler);
  }
  m_old_handler = NULL;
}

SigInstaller::~SigInstaller()
{
  release();
}

int gSignalSet = 0;
bool gSIGINT_is_fatal;

int pism_check_signal()
{
  int rv = gSignalSet;
  if (rv) {
    gSignalSet = 0;
  }
  return 0;
}

void pism_sigint_handler(int sig)
{
  if (sig == SIGINT) {
    if (gSIGINT_is_fatal) {
      throw pism::RuntimeError("caught signal SIGTERM.");
    } else {
      PetscPrintf(PETSC_COMM_WORLD, "\nCaught signal SIGTERM, waiting to exit.\n");
      gSignalSet = SIGINT;
    }
  }
}
