**********
invert.ssa
**********

Classes used for inverting the SSA from :math:`\tau_c` or hardness.  
The forward problem for inversion is defined by an :cpp:class:`IP_SSAForwardProblem` which knows about solving the SSA from a given value of :math:`\tau_c` or hardness. It needs to be initialized from a number of auxiliary data (geometry, grid, enthalpy converter, etc.)  The construction of all this data, along with the 
:cpp:class:`IP_SSATaucForwardProblem` is generally done by constructing 
a :class:`~PISM.invert.ssa.SSAForwardRunFromInputFile`, which is a subclass 
of :class:`PISM.ssa.SSARun`.

With the forward problem defined, a subclass of :class:`PISM.invert.ssa.InvSSASovler`
needs to be constructed that will apply a specific algorithm for solving 
the inverse problem.  The factory function :func:`~.createInvSSASolver`
constructs a solver based on command line arguments.

.. container:: custom-index

    .. raw:: html

        <script type="text/javascript" src='../../_static/pymunk.js'></script>

.. automodule:: PISM.invert.ssa
   :members:
   :undoc-members:
   :private-members:
   :special-members:
   :show-inheritance:
