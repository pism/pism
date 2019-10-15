.. include:: ../../global.txt

.. _sec-petscoptions:

PETSc options for PISM users
----------------------------

All PETSc programs including PISM accept command line options which control how PETSc
distributes jobs among parallel processors, how it solves linear systems, what additional
information it provides, and so on. The PETSc manual :cite:`petsc-user-ref` is the complete
reference on these options. We list some here that are useful to PISM users. They can be
mixed in any order with PISM options.

Both for PISM and PETSc options, there are ways of avoiding the inconvenience of long
commands with many runtime options. Obviously, and as illustrated by examples in the
previous sections, shell scripts can be set up to run PISM. But PETSc also provides two
mechanisms to give runtime options without retyping at each run command.

First, the environment variable ``PETSC_OPTIONS`` can be set. For example, a sequence of
runs might need the same refined grid, and you might want to know if other options are
read, ignored, or misspelled. Set (in Bash):

.. code-block:: none

   export PETSC_OPTIONS="-Mx 101 -My 101 -Mz 51 -options_left"

The runs

.. code-block:: none

   pismv -test F -y 100
   pismv -test G -y 100

then have the same refined grid in each run, and the runs report on which options were
read.

Alternatively, the file ``.petscrc`` is always read, if present, from the directory where
PISM (i.e. the PETSc program) is started. It can have a list of options, one per line. In
theory, these two PETSc mechanisms (``PETSC_OPTIONS`` and ``.petscrc``) can be used
together.

Now we address controls on how PETSc solves systems of linear equations, which uses the
PETSc "KSP" component (Krylov methods). Such linear solves are needed each time the
nonlinear SSA stress balance equations are used (e.g. with the option ``-stress_balance
ssa -ssa_method fd``).

Especially for solving the SSA equations with high resolution on multiple processors, it
is recommended that the option :opt:`-ssafd_ksp_rtol` be set lower than its default value
of :math:`10^{-5}`. For example,


.. code-block:: none

   mpiexec -n 8 ssa_testi -Mx 3 -My 769 -ssa_method fd

may fail to converge on a certain machine, but adding "``-ssafd_ksp_rtol 1e-10``" works
fine.

There is also the question of solver *type*, using option :opt:`-ssafd_ksp_type`. Based on
one processor evidence from ``ssa_testi``, the following are possible choices in the sense
that they work and allow convergence at some reasonable rate: ``cg``, ``bicg``, ``gmres``,
``bcgs``, ``cgs``, ``tfqmr``, ``tcqmr``, and ``cr``. It appears ``bicg``, ``gmres``,
``bcgs``, and ``tfqmr``, at least, are all among the best. The default is ``gmres``.

Actually the KSP uses preconditioning. This aspect of the solve is critical for parallel
scalability, but it gives results which are dependent on the number of processors. The
preconditioner type can be chosen with :opt:`-ssafd_pc_type`. Several choices are
possible, but for solving the ice stream and shelf equations we recommend only
``bjacobi``, ``ilu``, and ``asm``. Of these it is not currently clear which is fastest;
they are all about the same for ``ssa_testi`` with high tolerances (e.g. ``-ssafd_picard_rtol
1e-7`` ``-ssafd_ksp_rtol 1e-12``). The default (as set by PISM) is ``bjacobi``. To force
no preconditioning, which removes processor-number-dependence of results but may make the
solves fail, use ``-ssafd_pc_type none``.

For the full list of PETSc options controlling the SSAFD solver, run

.. code-block:: none

   ssa_testi -ssa_method fd -help | grep ssafd_ | less
