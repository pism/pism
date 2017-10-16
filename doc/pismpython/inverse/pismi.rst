.. _pismi:

``pismi.py``
==============


The ``pismi.py`` command-line tool is the primary interface for running SSA 
inversions.  It performs the following tasks:

1. Reads PISM model state variables and additional inverse data 
   from input files.

2. Optionally subtracts SIA velocities from observed surface velocities.

3. Runs a specified inversion algorithm.

4. Saves data and logs to an output NC file.


Data Files
----------

There are up to three NC files involved in the process:

* A model state file containing PISM's model state variables as generated
  by a spin-up run or otherwise.
* A inverse data file containing additional data needed for inversions
  (e.g. surface velocity measurements)
* An output file where the inverse solution (and other related variables) 
  are saved.

The model state file is specified with either :cfg:`-i` or :cfg:`-a`.  If
:cfg:`-a` is used, then ``pismi`` runs in append mode; the output file
is the same as the input file and new data is appended to the end of it
at the end of the solve.  Otherwise an output file must be specified with
:cfg:`-o` and the contents of the output file will be overwritten at the end
of the run.

The :cfg:`-inv_data` provides the name of the inverse data file. It can be the 
same as the input file.

The output file always contains enough data to run an inversion; ``pismi`` 
saves to it any variables from the model state and inverse data files
that were used in the course of the inversion.  Thus:

.. code-block:: none

   pismi.py -i model.nc -inv_data inv.nc -o solution.nc [OTHER FLAGS]
  
   pismi.py -i solution.nc -o solution2.nc [OTHER FLAGS]

should always result in a second inversion run that terminates immediately
because it is already at a solution.  The data in :file:`solution.nc` and
:file:`solution2.nc` should be identical.

.. _ObsSSAVel:

Observed SSA Velocities
-----------------------------

SSA velocities to be matched by the inversion algorithms are 
provided in the inverse data file, and can be provided directly
as the variables :ncvar:`u_ssa_observed`\ /:ncvar:`v_ssa_observed`.

On the other hand, the PISM hybrid model :cite:`BBssasliding` 
for horizontal surface velocities :math:`\vU` writes them as a sum

.. math::

  \vU = \vU_{\mathrm{SIA}} + \vU_{\mathrm{SSA}}.

where :math:`\vU_{\mathrm{SIA}}` and :math:`\vU_{\mathrm{SSA}}`
come from the SIA and SSA approximations.  Hence, if :math:`\vU_\obs`
are observed velocities at the surface, then an SSA inversion 
that is consistent with the PISM model should try to match

.. math::
  \vU_{\mathrm{SSA}} = \vU_\obs - \vU_{\mathrm{SIA}}.

If the inverse data file is missing 
:ncvar:`u_ssa_observed`\ /:ncvar:`v_ssa_observed` but contains
:ncvar:`u_surface_observed`\ /:ncvar:`v_surface_observed`, 
then ``pismi`` computes SIA velocities and 
subtracts them from the observed surface velocities
to arrive at :ncvar:`u_ssa_observed`\ /:ncvar:`v_ssa_observed`.

It may be the case that surface observations are not available at 
all grid points.  The variable :ncvar:`vel_misfit_weight` can be 
provided in the input file and can be used to indicate missing
values, or alternative weightings, as described in :ref:`state 
functionals <statefunc>`. If :ncvar:`vel_misfit_weight` is missing
it is assumed to be equal to 1 everywhere.

.. _pismi_design_var:

Design Variable
---------------

The inversion design variable is one of effective yield stress
:math:`\tau_c` or averaged hardness :math:`B`, and is specified
using :cfg:`-inv_ssa tauc` or :cfg:`-inv_ssa hardav` respectively.
The default is :cfg:`tauc`.

A parameterization for the design variable must also be specified
using :cfg:`-inv_design_param` and a name of one of the
:ref:`parameterizations <DesignParam>`.

The inversion algorithms require a best initial estimate for the
design variable, which is part of the *a-priori* data used to
regularize the inversion.  If there is a variable :ncvar:`tauc_prior`
in the inverse data file, it specifies the initial value.
Otherwise the initial value is taken from :ncvar:`tauc`
in the input file.  Use :cfg:`-no_use_tauc_prior` to
ignore the value of :ncvar:`tauc_prior` in the inverse data file
and to force the use of :ncvar:`tauc` instead.

Locations where :math:`\tau_c` (or :math:`B`) are to be held
constant at their initial estimates can be specified
with :ncvar:`zeta_fixed_mask`.  By default, these locations
are determined automatically.  If :ncvar:`zeta_fixed_mask`
is provided in the inverse data file, it will be used instead.
Use :cfg:`-no_use_zeta_fixed_mask` to disable the use
of :ncvar:`zeta_fixed_mask`.

In regions where PISM overrides the value of :math:`\tau_c` or
:math:`B` (i.e. in floating regions for :math:`\tau_c`) the 
initial estimate is adjusted to account for the PISM model.
**This might be a bad thing**.

At the end of inversion, the solution is saved as :ncvar:`tauc` in 
the output file.  Additionally, the final value of the
parameterized design variable :math:`\zeta` is saved as :ncvar:`zeta_inv`.

For hardness inversions, replace ``tauc`` with ``hardav``
in these variable and flag names.

Design and State Functionals
----------------------------

The choice of design and state functionals are made
with :cfg:`-inv_state_func` and :cfg:`-inv_design_func`
with a value among those documented in :ref:`state <statefunc>` and
:ref:`design <designfunc>` functional sections.


Inverse Algorithm Selection
---------------------------

The choice of inverse algorithm is made with the 
option :cfg:`-inv_method` with a value 
among those documented in the :ref:`iterative gradient <InvGradAlg>`
and :ref:`Tikhonov <TikhonovAlg>` algorithm sections.  The :cfg:`-inv_max_it`
flag determines the maximum number of iterations allowed by the algorithm.


Regularization Constants
------------------------

For iterative gradient algorithms,
:cfg:`-inv_target_misfit` specifies
the :ref:`stopping criterion <InvGradStop>`.

For Tikhonov algorithms use :cfg:`-tikhonov_penalty`
to specify the :ref:`penalty parameter <TikhonovAlg>`.

See also the discussion on :ref:`Tikhonov minimization
convergence <TikConverge>`.

Other SSA-Related Flags
-----------------------

Any flags that affect the SSA in a usual PISM run need to be
specified for ``pismi`` as well.  These include, but are not
limited to,

* :cfg:`-ssa_dirichlet_bc`\ : Apply Dirichlet boundary conditions.
* :cfg:`-regional`\ : Use PISM regional model semantics.
* :cfg:`-pseudo_plastic`\ : Use the pseudo-plastic till model.
* :cfg:`-pseudo_plastic_q`\ : Sets the value of :math:`q` for the pseduo plastic till model.
* :cfg:`-flow_law`\ : Sets the ice flow law model (e.g. Patterson-Budd polythermal Glen ice via :cfg:`pb`).

Model State File Contents
-------------------------

The model state file must contain the following variables:

  1. Bedrock elevation :ncvar:`topg`
  2. Ice thickness :ncvar:`thk`
  3. Enthalpy :ncvar:`enthalpy`
  
If Dirichlet boundary conditions are being used (:cfg:`-ssa_dirichlet_bc`),
the model state file must contain

  4. SSA Dirichlet velocities :ncvar:`vel_ssa_bc`
  5. Dirichlet mask :ncvar:`bc_mask` specifying where Dirichlet conditions 
     apply.
  
If PISM is being used in regional model mode (:cfg:`-regional`), this last variable is replaced with
  
  5. :ncvar:`no_model_mask`.

An initial estimate for the design variable :ncvar:`tauc` or :ncvar:`hardav`
can be provided as well, as discussed in the :ref:`design variable <pismi_design_var>` section.

Inverse Data File Contents
--------------------------

The following variables may be present in the inverse data file:

  1. :ncvar:`u_ssa_observed`\ /:ncvar:`v_ssa_observed`: Target SSA
     velocities to be matched by the inversion algorithm.

  2. :ncvar:`u_surface_observed`\ /:ncvar:`v_surface_observed`: Observed
     surface velocities used to generate 
     :ncvar:`u_ssa_observed`\ /:ncvar:`v_ssa_observed`: 
     :ref:`if needed <ObsSSAVel>`.

  3. :ncvar:`vel_misfit_weight`\ : The weight function discussed
     in the :ref:`oberved SSA Velocity <ObsSSAVel>` section.

  4. :ncvar:`tauc_prior` or :ncvar:`hardav_prior`\ : The
     *a-priori* best estimate for the physical design variable, overriding the
     value in the model state file.  

  5. :ncvar:`zeta_fixed_mask`\ : Locations where the design variable is 
      to be  held at its initial estimate.  If this variable is not present,
      an appropriate mask will be generated, unless
      :cfg:`-no_use_zeta_fixed_mask` is specified.

  6. :ncvar:`zeta_inv`\ : The initial value of the parameterized 
     design variable to start iterating from.  If it is absent,
     it will be constructed from :ncvar:`tauc_prior` 
     (or :ncvar:`hardav_prior`).

All of these are optional, except:

  * At least one of :cfg:`ssa_observed` or :cfg:`surface_observed`
    velocities must be present, with :cfg:`ssa_observed` velocities
    used if both are present.

  * For :math:`\tau_c` inversions, if :ncvar:`tauc_prior` 
    is not present, or if :cfg:`-no_inv_use_tauc_prior` is set,
    then :ncvar:`tauc` must be present in the input file.  A similar
    caveat holds for hardness inversions replacing :ncvar:`tauc`.
    with :ncvar:`hardav`.


Output File
-----------

The following variables are written to the output file,
in addition to a number of variables that were provided
in the model state and inverse data files:

  * :ncvar:`tauc` or :ncvar:`hardav`\ : The value of
    the design variable solved for by inversion.
  * :ncvar:`zeta_inv`: The last computed value of the
    parameterized design variable :math:`\zeta`.
  * :ncvar:`u_ssa_inv`\ /:ncvar:`v_ssa_inv` : The
    SSA velocities corresponding to the design
    variable arrived at by inversion.
  * :ncvar:`u_inv_ssa_residual`\ :ncvar:`v_inv_ssa_residual`:
    The difference between observed SSA velocities and the
    velocities arrived at by inversion.
  * :ncvar:`inv_ssa_residual` : The magnitude of the velocity
    residuals.

The output file also contains a log of the inversion run
in the NC variable :ncvar:`pismi_log`.  The output file 
also contains a log of the misfit at each iteration
in the variable :ncvar:`inv_ssa_misfit`.  For 
:cfg:`-inv_state_func meansquare`, the values will
be square roots of the misfit functional, in units of m/a.  
Otherwise, these will be the values of the misfit functional
itself.

A copy of the command line used to run the inversion is saved
in the :ncvar:`history` attribute of the output file.

Prep File and Listeners
-----------------------

A python module can be provided to perform additional setup
prior to starting inversion.  Use the :cfg:`-inv_prep_module`
to indicate a python module containing a function
:func:`prep_solver(solver)`, which receives 
a :class:`PISM.invert.ssa.Solver` object as its argument.
To attach a listener object to be called at each iteration,
use :func:`solver.addIterationListener`.

See also the :ref:`listener <Listeners>` documentation.


Restarting Inversion
--------------------

At each iteration of the inversion, a copy of the current
parameterized design variable :math:`\zeta` is saved
as :ncvar:`zeta_inv` in the output file.  If for some reason
``pismi`` is interrupted (e.g. control-C), inversion can
be restarted from the last saved iterate by specifying
:cfg:`-inv_restart` along with all of the other 
command-line flags used originally to run the inversion.

