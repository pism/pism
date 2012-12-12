#ifndef _PISM_SIGNAL_H
#define _PISM_SIGNAL_H

int pism_check_signal();
void pism_sigint_handler(int sig);

extern bool gSIGINT_is_fatal;

//! Installs a signal handler on construction; deinstalls on destruction.
class SigInstaller
{
public:
  //! Installs handle \a new_handler for signal \a sig.
  SigInstaller(int sig, void (*new_handler)(int));
  //! Restores the signal handler to its previous value.
  void release();
  //! 
  ~SigInstaller();
private:
  void (*m_old_handler)(int);
  int m_sig;
};

#endif
