// Copyright (C) 2011 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "PSExternal.hh"

// Note: EBM driver code should NOT use any PETSc calls! (The process running
// this code does not initialize PETSc.)

EBM_driver::EBM_driver(MPI_Comm my_com) {
  inter_comm = my_com;
  command = "WILL BE INITIALIZED LATER";
}

int EBM_driver::run() {
  int ierr;
  bool ready;
  int done = 0,
    wait_counter = 0,
    wait_message_counter = 1;

  // Initialization: get the command.
  {
    char tmp[PETSC_MAX_PATH_LEN];
    MPI_Recv(tmp, PETSC_MAX_PATH_LEN, MPI_CHAR, 0, TAG_EBM_COMMAND, inter_comm, NULL);
    command = tmp;
  }

  while ( !done ) {
    char filename[PETSC_MAX_PATH_LEN];

    MPI_Status status;
    int flag;
    MPI_Iprobe(0, MPI_ANY_TAG, inter_comm, &flag, &status);

    if ( !flag ) {
      // no message
      double sleep_interval = 0.01,
        threshold = 3600,
        message_interval = 5;  // print a message every 5 seconds
      struct timespec rq;
      rq.tv_sec = 0;
      rq.tv_nsec = (long)(sleep_interval*1e9); // convert to nanoseconds

      if (wait_counter * sleep_interval > threshold) {
        fprintf(stderr, "ERROR: spent %1.1f minutes waiting for PISM... Giving up...\n",
                threshold / 60.0);

        done = 1;
        continue;
      }

      // waiting...
      if (sleep_interval * wait_counter / message_interval  > wait_message_counter) { // 
        fprintf(stderr, "EBM driver: Waiting for a message from PISM... (%1.1f seconds so far)\n",
                message_interval * wait_message_counter);
        wait_message_counter++;
      }
      nanosleep(&rq, 0);
      wait_counter++;
      continue;
    }

    // we got a message; reset counters
    wait_counter = 0;
    wait_message_counter = 1;

    switch (status.MPI_TAG) {
    case TAG_EBM_RUN:
      {
        int tmp, ebm_status;
        MPI_Recv(&tmp, 1, MPI_INT, 0, TAG_EBM_RUN, inter_comm, NULL);
        
        ierr = run_ebm();

        if (ierr != 0) {
          // the run failed
          ebm_status = EBM_STATUS_FAILED;
          done = 1;
        } else {
          // we're OK
          ebm_status = EBM_STATUS_READY;
        }

        MPI_Send(&ebm_status, 1, MPI_INT, 0, TAG_EBM_STATUS, inter_comm);

        break;
      }
    case TAG_EBM_STOP:
      {
        int tmp;
        MPI_Recv(&tmp, 1, MPI_INT, 0, TAG_EBM_STOP, inter_comm, NULL);
        done = 1;
        break;
      }
    default:
      {
        fprintf(stderr, "ERROR: Unknown message tag: %d\n", status.MPI_TAG);
        break;
      }
    }
  } // while (not done)

  return 0;
}

int EBM_driver::run_ebm() {
  int ierr;

  fprintf(stderr, "EBM driver: running '%s'...\n", command.c_str());

  ierr = system(command.c_str());

  if (ierr != 0) {
    fprintf(stderr, "EBM driver: command exited with the return code %d\n", ierr);
    return ierr;
  } else {
    fprintf(stderr, "EBM driver: command executed succesfully.\n", ierr);
  }

  return 0; 
}

