.. include:: ../../global.txt

Initialization from a saved model state
---------------------------------------

"Initialization" has the specific, simple meaning in PISM that option ":opt:`-i`" was
used. If a previous PISM run has saved a NetCDF file using ":opt:`-o`" then that file will
contain complete initial conditions for continuing the run. The output file from the last
run can be loaded with ":opt:`-i`":

.. code-block:: none

   pisms -eisII A -y 100 -o foo.nc
   pisms -eisII A -i foo.nc -y 100 -o bar.nc

As noted verification tests (section :ref:`sec-verif`) and simplified-geometry experiments
(section :ref:`sec-simp`) do not need input files at all because they initialize from
formulas in the source code. They can, however, be continued from saved model states using
:opt:`-i`. Specifying the simplified geometry experiment or verification test *is*,
however, necessary if the run is to continue with the climate inputs for that experiment
or test. For example, based on the above ``pisms`` runs, it is valid to do

.. code-block:: none

   pismr -i foo.nc -y 100 -o bar.nc

but the climate and other parameters use PISM default values, and thus are not
(necessarily) the values specified in EISMINT II.

As a technical note about saved states, a PISM run with :opt:`-stress_balance ssa` also
saves the last SSA velocities to the output file in variables :var:`u_ssa` and
:var:`v_ssa`. The presence of these velocities adds efficiency in restarting because an
initial estimate speeds up the solution of the SSA stress balance equations. If you want
to use :opt:`-i` but also ignore these velocities then set
:config:`stress_balance.ssa.read_initial_guess` to *false*.

.. _sec-i-format:

``-i`` file format
^^^^^^^^^^^^^^^^^^

PISM produces CF-1.5 compliant NetCDF files. The easiest way to learn the output format
*and* the :opt:`-i` format is to do a simple run and then look at the metadata in the
resulting file, like this:

.. code-block:: none

   pisms -eisII A -y 10 -o foo.nc
   ncdump -h foo.nc | less

Note that variables in the output file have a ``pism_intent`` attribute. When
``pism_intent`` is ``diagnostic``, the variable can be deleted from the file without
affecting whether PISM can use it as a :opt:`-i` input file. Variables with
``pism_intent`` is ``model_state``, by contrast, must be present when using :opt:`-i`.

The automatically-produced :var:`time` variable has a ``units`` attribute like ``"seconds
since 1-1-1"`` because the CF metadata conventions require a reference date.

.. FIXME: double-check the statement below

By default PISM ignores this reference date except when it is used in unit conversions
based on a calendar (see below).
