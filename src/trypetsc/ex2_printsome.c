
/* Program usage:  mpiexec ex2 [-help] [all PETSc options] */

static char help[] = "Synchronized printing.\n"
                     "Modification by Bueler, of src/sys/examples/tutorials/ex2.c,\n"
                     "which shows how to get output even if one of the processors\n"
                     "in the group doesn't call the print statement at all.\n\n";

/*T
   Concepts: introduction to PETSc;
   Concepts: printing^synchronized
   Concepts: printing^in parallel
   Concepts: printf^synchronized
   Concepts: printf^in parallel
   Processors: n
T*/
 
#include "petsc.h"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,size;

  /*
    Every PETSc program should begin with the PetscInitialize() routine.
    argc, argv - These command line arguments are taken to extract the options
                 supplied to PETSc and options supplied to MPI.
    help       - When PETSc executable is invoked with the option -help, 
                 it prints the various options that can be applied at 
                 runtime.  The user can use the "help" variable place
                 additional help messages in this printout.
  */
  ierr = PetscInitialize(&argc,&argv,PETSC_NULL,help);CHKERRQ(ierr);

  /* 
     The following MPI calls return the number of processes
     being used and the rank of this process in the group.
   */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  /* 
     Here we would like to print only one message that represents
     all the processes in the group.  We use PetscPrintf() with the 
     communicator PETSC_COMM_WORLD.  Thus, only one message is
     printed representing PETSC_COMM_WORLD, i.e., all the processors.
  */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of processors = %d, rank = %d\n",size,rank);CHKERRQ(ierr);

  /* 
     Here we would like to print info from each process, such that
     output from process "n" appears after output from process "n-1".
     To do this we use a combination of PetscSynchronizedPrintf() and
     PetscSynchronizedFlush() with the communicator PETSC_COMM_WORLD.
     All the processes print the message, one after another. 
     PetscSynchronizedFlush() indicates that the current process in the
     given communicator has concluded printing, so that the next process
     in the communicator can begin printing to the screen.
  */
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] Synchronized Hello World.\n",rank);CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] Synchronized Hello World - Part II.\n",rank);CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);

  /*
    Here a barrier is used to separate the two stages.
  */
  ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);

  /*
    Here we simply use PetscPrintf() with the communicator PETSC_COMM_SELF
    (where each process is considered separately).  Thus, this time the
    output from different processes does not appear in any particular order.
    
    Unlike in ex2.c we show that one process can just not write at all and
    still we will get output.
  */
  if (rank != 1)   {
    ierr = PetscPrintf(PETSC_COMM_SELF,
       "[%d] Hello World from PETSC_COMM_SELF (except rank==1)\n",rank);CHKERRQ(ierr);
  }

  /*
    Here a barrier is used to separate the two states.
  */
  ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);


  /*
    This is also o.k. but the synchronization is effectively meaningless.
    At least it does not cause hang.  This would cause a hang if we used 
    "PetscSynchronizedPrintf(PETSC_COMM_WORLD,..." instead.
  */
  if (rank != 1)   {
    ierr = PetscSynchronizedPrintf(PETSC_COMM_SELF,
       "[%d] Hello World from PETSC_COMM_SELF meaninglessly synchronized (except rank==1)\n",rank);CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_SELF);CHKERRQ(ierr);
  }

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).  
     See the PetscFinalize() manpage for more information.
  */
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
