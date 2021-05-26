.. include:: ../../global.txt

.. _sec-flowline-modeling:

Using PISM for flow-line modeling
---------------------------------

As described in sections :ref:`sec-coords` and :ref:`sec-grid`, PISM is a
three-dimensional model. Moreover, parameters :config:`grid.Mx` and :config:`grid.My` have
to be greater than or equal to three, so it is not possible to turn PISM into a 2D
(flow-line) model by setting :config:`grid.Mx` or :config:`grid.My` to `1`.

There is a way around this, though: by using :config:`grid.periodicity` to tell
PISM to make the computational grid `y`\-periodic and providing initial and boundary
conditions that are functions of `x` only one can ensure that there is no flow in
the `y`\-direction.

In this case :config:`grid.My` can be any number; we want to avoid unnecessary
computations, though, so :config:`grid.My` of `3` is the obvious choice.

One remaining problem is that PISM still expects input files to contain both ``x`` and
``y`` dimensions. To help with this, PISM comes with a Python script ``flowline.py`` that
turns NetCDF files with `N` grid points along a flow line into files with 2D fields
containing `N\times3` grid points.\ [#]_

Here's an example which uses the script ``examples/preprocessing/flowlineslab.py`` to
create a minimal, and obviously unrealistic, dataset. The file ``slab.nc`` created by this
script contains all the required information but is not ready to use with PISM. Proceed as
follows:

.. code-block:: bash

   examples/preprocessing/flowlineslab.py # creates slab.nc with only an x-direction
   flowline.py -o slab-in.nc --expand -d y slab.nc

This produces a PISM-ready ``slab-in.nc``. Specifically, ``flowline.py`` "expands" its
input file in the `y`\-direction. Now we can "bootstrap" from ``slab-in.nc``:

.. code-block:: none

   mpiexec -n 2 pismr \
           -surface given \
           -bootstrap -i slab-in.nc \
           -Mx 201 -Lx 1000 \
           -My 3 -Ly 4 -periodicity y \
           -Lz 2000 -Mz 11 \
           -y 10000 -o pism-out.nc

To make it easier to visualize data in the file created by PISM, "collapse" it using NCO_:

.. code-block:: none

   ncks -O -d y,1 pism-out.nc slab-out.nc
   ncwa -O -a time,y slab-out.nc slab-out.nc

.. rubric:: Footnotes

.. [#] This script requires the ``numpy`` and ``netCDF4`` Python modules. Run
       ``flowline.py --help`` for a full list of options.
