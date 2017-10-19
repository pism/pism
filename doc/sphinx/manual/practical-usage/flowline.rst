.. include:: ../../global.txt

.. _sec-flowline-modeling:

Using PISM for flow-line modeling
---------------------------------

As described in sections :ref:`sec-coords` and :ref:`sec-grid`, PISM is a
three-dimensional model. Moreover, parameters ``Mx`` and ``My`` have to be greater than or
equal to three, so it is not possible to turn PISM into a 2D (flow-line) model by setting
``Mx`` or ``My`` to 1.

There is a way around this, though: by using the :opt:`-periodicity` option to tell PISM
to make the computational grid :math:`y`-periodic and providing initial and boundary
conditions that are functions of :math:`x` only one can ensure that there is no flow in
the :math:`y`\-direction. (Option :opt:`-periodicity` takes an argument specifying the
direction: ``none``, ``x``, ``y`` and ``xy`` --- for "periodic in both X- and
Y-directions".)

In this case ``Mx`` can be any number; we want to avoid unnecessary computations, though,
so "``-Mx 3``" is the obvious choice.

One remaining problem is that PISM still expects input files to contain both ``x`` and
``y`` dimensions. To help with this, PISM comes with a Python script ``flowline.py`` that
turns NetCDF files with :math:`N` grid points along a flow line into files with 2D fields
containing :math:`N\times3` grid points.\ [#]_

Here's an example which uses the script ``util/flowlineslab.py`` to create a minimal, and
obviously unrealistic, dataset. A file ``slab.nc`` is created by ``util/flowlineslab.py``,
but it is not ready to use with PISM. Proceed as follows, after checking that ``util/`` is
on your path:

.. code-block:: bash

   flowlineslab.py                         # creates slab.nc with only an x-direction
   flowline.py -o slab-in.nc --expand -d y slab.nc


produces a PISM-ready ``slab-in.nc``. Specifically, ``flowline.py`` "expands" its input
file in the y-direction. Now we can "bootstrap" from ``slab-in.nc``:

.. code-block:: none

   mpiexec -n 2 pismr -surface given -i slab-in.nc -bootstrap -periodicity y \
           -Mx 201 -My 3 -Lx 1000 -Ly 4 -Lz 2000 -Mz 11 -y 10000 -o pism-out.nc

To make it easier to visualize data in the file created by PISM, "collapse" it:

.. code-block:: none

   flowline.py -o slab-out.nc --collapse -d y pism-out.nc

.. rubric:: Footnotes

.. [#] This script requires the ``numpy`` and ``netCDF4`` Python modules. Run
       ``flowline.py --help`` for a full list of options.
