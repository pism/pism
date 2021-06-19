=====================
Inverse Code Overview
=====================

This document is an brief
overview of how the algorithms for solving inverse problems have been grafted onto the main PISM library, intended to
help with future maintenance of the code.

Inverse problem solving in PISM consists of the following components:

  * ``siple`` An older python-based library, that includes iterative 
    gradient algorithms for solutions of inverse problems, including 
    steepest descent, nonlinear conjugate gradients, and the so-called 
    incomplete Gauss-Newton method.

  * A newer C++ codebase within PISM implementing Tikhonov-type 
    regularization algorithms, many of which use the
    TAO :cite:`tao-user-ref` optimization library.

  * A python PISM library that manages the interface between standard 
    PISM components, selection of inverse algorithms, and the two inverse
    libraries mentioned above.

A typical inversion run for reconstructing basal yield stresses from surface velocities starts with the script :file:`pismi.py`. It instantiates an instance of :class:`PISM.inv_ssa.InvSSARun`, which knows all of the data required to compute solutions of the SSA, and knows about the map from basal yield stress to surface velocities, but knows nothing more about the details 
of the inverse problem, how it will be regularized, or the algorithm that will
be used to solved the regularized problem.  The :class:`~PISM.inv_ssa.InvSSARun` instance is then handed by 
:file:`pismi.py` to a solver, e.g. 
:class:`PISM.inv_ssa_tao.InvSSASolver_Tikhonov` or 
:class:`PISM.inv_ssa_siple.InvSSASolver_Iterative`.  The solvers are responsible
for wiring up the algorithm library specific components to the :class:`PISM.inv_ssa.InvSSARun`, setting up any details about the regularization or algorithm, and present a uniform interface to call to solve the inverse problem. The actual solution of the problem is then delegated by the solver to the underlying library.

Here are the key players common to all inversions.

  * :file:`pismi.py`\ : user level script 
  * :cpp:class:`SSAFEM`\ : Base class for solving the SSA via 
    finite elements using PETSc SNES.
  * :class:`PISM.model.ModelData`\ : Contains vectors and physics shared 
    between SSA and SIA solvers.
  * :class:`PISM.ssa.SSARun`\ : Responsible for putting together
    a :class:`~PISM.model.ModelData` and a 
    :cpp:class:`SSAFEM` (or related class)  for solving the SSA.
    It can then be used to generate solutions of the SSA.
  * :class:`PISM.sia.computeSIASurfaceVelocities`\ : Utility function
    for solving the SIA from a :class:`~PISM.model.ModelData`.
  * :cpp:class:`IP_SSATaucForwardProblem`\ : Subclass of
    :cpp:class:`SSAFEM` that contains additional methods for 
    defining a forward problem (i.e. a map from :math:`\tau_c` to 
    SSA velocities).  In particular, it contains code for various
    related linearizations.  
  * :class:`PISM.inv_ssa.SSAInvRun`\ : Subclass of :class:`~PISM.ssa.SSARun`
    that reads in additional inverse data and 
    where the underlying :cpp:class:`SSAFEM` 
    is a :cpp:class:`IP_SSATaucForwardProblem`.

For TAO-based Tikhonov problems, we have the additional parties involved.

  * :cpp:class:`TaoBasicSolver`\ : C++ interface to TAO optimization routines.
  * :cpp:class:`IPTaoTikhonovProblem`\ : Defines a generic Tikhonov problem to
    be solved by a :cpp:class:`TaoBasicSolver`.
  * :cpp:class:`IP_SSATaucTaoTikhonovProblem`\ : Subclass of 
    :cpp:class:`IPTaoTikhonovProblem` for SSA Tikhonov problems. Uses a
    :cpp:class:`IP_SSATaucForwardProblem` to define the forward problem.
  * :class:`PISM.inv_ssa_tao.InvSSASolver_Tikhonov`\ : Python glue
    between :cpp:class:`IPTaoTikhonovProblem` and ``pismi``.

For ``siple``-based inversions, the following classes are important.

  * :class:`siple.gradient.forward.NonlinearForwardProblem`\ :
    Defines a ``siple`` forward problem wrapping a :cpp:class:`IP_SSATaucForwardProblem`.
  * :class:`siple.gradient.InvertNLCG`/:class:`siple.gradient.InvertIGN`\ :
    Base classes for iterative gradient solvers.
  * :class:`PISM.inv_ssa_siple.InvertSSANLCG`/:class:`PISM.inv_ssa_siple.InvertSSAIGN`\ : SSA-based iterative gradient solvers.  
  * :class:`PISM.inv_ssa_siple.InvSSASolver_Gradient`\ : Python glue
    between ``siple`` and ``pismi``.
