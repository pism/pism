.. _Listeners:

Iteration Listeners
===================

The PISM inverse code provides a flexible mechanism for providing user
callbacks that can report diagnostic information or perform other activities
at each iteration of an inversion algorithm.

A listener is either a function with signature::

  def my_listener( solver, iteration, data ):

or a callable class that implements::

  def __call__(self, solver, iteration, data):

where

  * ``solver`` is an instance of :class:`PISM.inverse.InvSSASolver`
  * ``iteration`` is the current iteration count
  * ``data`` is a :class:`PISM.util.Bunch` with data as described below.

Prior to starting the solution algorithm, listeners are attached to the 
solver using :meth:``PISM.inverse.InvSSASolver.addIterationListener``.
Each listener is called in turn at each step of the inversion,
with extra data provided in the ``data`` argument.  The following
attributes are always available in ``data``:

  * ``zeta``: The current value of the parameterized design variable.
  * ``u``: The current state value.
  * ``zeta_step``: The most recent adjustment to ``zeta``.
  * ``residual``: The difference between the state value ``u`` and the
    target state value.

Tikhonov algorithms may also provide:

  * ``JDesign``: The value of the design functional.
  * ``JState``: The value of the state functional.
  * ``grad_JDesign``: The derivative of the design functional
                      with respect to the (parameterized) design variable.
  * ``grad_JState``: The derivative of the state functional
    with respect to the (parameterized) design variable.
  * ``grad_JTikhonov``: The derivative of the combined Tikhonov functional
    with respect to the (parameterized) design variable.
  * ``tikhonov_penalty``: The Tikhonov penalty parameter :math:`\eta`.

Iterative gradient algorithms may also provide:

  * ``target_misfit``: The desired value of the state misfit 
    :math:`\sqrt{J_S}`\ .
  * ``T_zeta_step``: The linearization of the forward map applied 
    to ``zeta_step``.
  * ``TStar_residual``: The adjoint of the linearization of the forward 
    map applied to ``residual``.

One can verify a desired variable is present and access it with the 
following idiom::

  JDesign = None
  if data.has_key('JDesign'):
    JDesign = data.JDesign

Vectors distributed across processors can be gathered to processor zero using
:class:`PISM.vec.ToProcZero`, which can be helpful for, e.g., plotting them.

.. _customListener:

Custom Listeners
----------------

To add a listener to an inversion run using ``pismi.py``, use
:cfg:`-inv_prep_module` to indicate a python module containing 
a function :func:`prep_solver(solver)`, which receives 
a :class:`PISM.invert.ssa.Solver` object as its argument.
Then call :func:`PISM.invert.ssa.Solver.addIterationListener`
as desired to attach the custom listeners. E.g. ::

  import PISM
  def myListener(solver, iteration, data):
    if data.has_key('JDesign'):
      PISM.logging.logMessage("At iteration %d, design functional value is: %g\n" 
        % (iteration, data.JDesign) )
    else:
      # Ignore
      pass

  
  def prep_solver(solver):
    solver.addIterationListener(myListener)
