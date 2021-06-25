.. include:: ../../global.txt

.. _sec-regridding:

Regridding
----------

.. FIXME: Mention that -regrid_file without -regrid_vars re-grids model state variables
   for all selected sub-models.

It is common to want to interpolate a coarse grid model state onto a finer grid or vice
versa. For example, one might want to do the EISMINT II experiment on the default grid,
producing output ``foo.nc``, but then interpolate both the ice thickness and the
temperature onto a finer grid. The basic idea of "regridding" in PISM is that one starts
over from the beginning on the finer grid, but one extracts the desired variables stored
in the coarse grid file and interpolates these onto the finer grid before proceeding with
the actual computation.

The transfer from grid to grid is reasonably general --- one can go from coarse to fine or
vice versa in each dimension `x,y,z` --- but the transfer must always be done by
*interpolation* and never *extrapolation*. (An attempt to do the latter will always
produce a PISM error.)

Such "regridding" is done using the :opt:`-regrid_file` and :opt:`-regrid_vars` commands
as in this example: }

.. code-block:: none

    pismr -eisII A -Mx 101 -My 101 -Mz 201 -y 1000 \
          -regrid_file foo.nc -regrid_vars thk,temp -o bar.nc

By specifying regridded variables "``thk,temp``", the ice thickness and temperature values
from the old grid are interpolated onto the new grid. Here one doesn't need to regrid the
bed elevation, which is set identically zero as part of the EISMINT II experiment A
description, nor the ice surface elevation, which is computed as the bed elevation plus
the ice thickness at each time step anyway.

A slightly different use of regridding occurs when "bootstrapping", as described in
section :ref:`sec-initboot` and illustrated by example in section :ref:`sec-start`.

See :numref:`tab-regridvar` for the regriddable variables using ``-regrid_file``.
Only model state variables are regriddable, while climate and boundary data generally are
not explicitly regriddable. (Bootstrapping, however, allows the same general interpolation
as this explicit regrid.)

.. list-table:: Regriddable variables.  Use ``-regrid_vars`` with these names.
   :header-rows: 1
   :name: tab-regridvar

   * - Name
     - Description
   * - :var:`age`
     - age of ice
   * - :var:`bwat`
     - effective thickness of subglacial melt water
   * - :var:`bmelt`
     - basal melt rate
   * - :var:`dbdt`
     - bedrock uplift rate
   * - :var:`litho_temp`
     - lithosphere (bedrock) temperature
   * - :var:`mask`
     - grounded/dragging/floating integer mask, see section :ref:`sec-floatmask`
   * - :var:`temp`
     - ice temperature
   * - :var:`thk`
     - land ice thickness
   * - :var:`topg`
     - bedrock surface elevation
   * - :var:`enthalpy`
     - ice enthalpy

Here is another example: suppose you have an output of a PISM run on a fairly coarse grid
(stored in ``foo.nc``) and you want to continue this run on a finer grid. This can be done
using ``-regrid_file`` along with ``-bootstrap``:

.. code-block:: none

   pismr -i foo.nc -bootstrap -Mx 201 -My 201 -Mz 21 -Lz 4000 \
         -regrid_file foo.nc -regrid_vars litho_temp,enthalpy -y 100 -o bar.nc \
         -surface constant

In this case all the model-state 2D variables present in ``foo.nc`` will be interpolated
onto the new grid during bootstrapping, which happens first, while three-dimensional
variables are filled using heuristics mentioned in section :ref:`sec-initboot`. Then
temperature in bedrock (``litho_temp``) and ice enthalpy (``enthalpy``) will be
interpolated from ``foo.nc`` onto the new grid during the regridding stage, overriding
values set at the bootstrapping stage. All of this, bootstrapping and regridding, occurs
before the first time step.

By default PISM checks the grid overlap and stops if the current computational domain is
not a subset of the one in a ``-regrid_file``. It is possible to disable this check and
allow constant extrapolation: use the option :opt:`-allow_extrapolation`.

For example, in a PISM run the ice thickness has to be lower than the vertical extent of
the computational domain. If the ice thickness exceeds ``Lz`` PISM saves the model state
and stops with an error message.

.. code-block:: none

   pismr -i input.nc -bootstrap -Mz 11 -Lz 1000 -z_spacing equal \
         -y 3e3 \
         -o too-short.nc
   PISM ERROR: Ice thickness exceeds the height of the computational box (1000.0000 m).
               The model state was saved to 'too-short_max_thickness.nc'.
               To continue this simulation, run with
               -i too-short_max_thickness.nc -bootstrap -regrid_file too-short_max_thickness.nc \
               -allow_extrapolation -Lz N [other options]
               where N > 1000.0000.

Regridding with extrapolation makes it possible to extend the vertical grid and continue a
simulation like this one --- just follow the instructions provided in the error message.
