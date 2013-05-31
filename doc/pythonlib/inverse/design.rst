===============================================
Overview of inversion components in PISM
===============================================

This document is an overview of how the algorithms for solving inverse problems have been grafted onto the main PISM library, intended to
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

A typical inversion run for reconstructing basal yield stresses from surface velocities starts with the script :file:`vel2tauc.py`. It instantiates an instance of :class:`PISM.inv_ssa.InvSSARun`, which knows all of the data required to compute solutions of the SSA, and knows about the map from basal yield stress to surface velocities, but knows nothing more about the details 
of the inverse problem, how it will be regularized, or the algorithm that will
be used to solved the regularized problem.  The :class:`~PISM.inv_ssa.InvSSARun` instance is then handed by :file:`vel2tauc.py` to a solver, e.g. 
:class:`PISM.inv_ssa_tao.InvSSASolver_Tikhonov` or 
:class:`PISM.inv_ssa_siple.InvSSASolver_Iterative`.  The solvers are responsible
for wiring up the algorithm library specific components to the :class:`PISM.inv_ssa.InvSSARun`, setting up any details about the regularization or algorithm, and present a uniform interface to call to solve the inverse problem. The actual solution of the problem is then delegated by the solver to the underlying library.

Here are the key players common to all inversions.(TODO: All the details!)

  * :file:`pismi.py`
  * :cpp:class:`SSAFEM`
  * :class:`PISM.model.ModelData`
  * :class:`PISM.ssa.SSARun`
  * :class:`PISM.sia.computeSIASurfaceVelocities`
  * :class:`PISM.inv_ssa.SSAInvRun`
  * :cpp:class:`IP_SSATaucForwardProblem`

For TAO-based Tikhonov problems, we have the additional parties involved.

  * :class:`PISM.inv_ssa_tao.InvSSASolver_Tikhonov`
  * :cpp:class:`TaoBasicSolver`
  * :cpp:class:`IPTaoTikhonovProblem`
  * :cpp:class:`IP_SSATaucTaoTikhonovProblem`

For ``siple``-based inversions, the following classes are important.

  * :class:`PISM.inv_ssa_siple.InvSSASolver_Iterative`
  * :class:`siple.gradient.forward.NonlinearForwardProblem`
  * :class:`siple.gradient.InvertNLCG`/:class:`siple.gradient.InvertIGN`
  * :class:`PISM.inv_ssa_siple.InvertSSANLCG`/:class:`PISM.inv_ssa_siple.InvertSSAIGN`

SOMETHING ABOUT LISTENERS
