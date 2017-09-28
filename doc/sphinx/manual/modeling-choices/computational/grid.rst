.. include:: ../../../global.rst

.. FIXME: add a section about the grid extent and the way to specify it (cell centers vs.
   cell corners, etc)

.. _sec-grid:

Spatial grid
------------

The PISM grid covering the computational box is equally spaced in horizontal (`x`
and `y`) directions. Vertical spacing in the ice is quadratic by default but
optionally equal spacing can be chosen; choose with options :opt:`-z_spacing`
[``quadratic``, ``equal``\] at bootstrapping. The grid read from a "``-i``" input file is
used as is. The bedrock thermal layer model always uses equal vertical spacing.

The grid is described by four numbers, namely the number of grid points ``Mx`` in the
`x` direction, the number ``My`` in the `y` direction, the number ``Mz`` in
the `z` direction within the ice, and the number ``Mbz`` in the `z` direction
within the bedrock thermal layer. These are specified by options :opt:`-Mx`, :opt:`-My`,
:opt:`-Mz`, and :opt:`-Mbz`, respectively. The defaults are 61, 61, 31, and 1,
respectively. Note that ``Mx``, ``My``, ``Mz``, and ``Mbz`` all indicate the number of
grid *points* so the number of grid *spaces* are one less, thus 60, 60, 30, and 0 in the
default case.

The lowest grid point in a column of ice, at `z=0`, coincides with the highest grid point
in the bedrock, so ``Mbz`` must always be at least one. Choosing ``Mbz`` `>1` is required
to use the bedrock thermal model. When a thermal bedrock layer is used, the distance
``Lbz`` is controlled by the :opt:`-Lbz` option. Note that ``Mbz`` is unrelated to the bed
deformation model (glacial isostasy model); see section :ref:`sec-beddef`.

In the quadratically-spaced case the vertical spacing near the ice/bedrock interface is
about four times finer than it would be with equal spacing for the same value of ``Mz``,
while the spacing near the top of the computational box is correspondingly coarser. For a
detailed description of the spacing of the grid, see the documentation on
``IceGrid::compute_vertical_levels()`` in the `PISM class browser <pism-browser_>`_.

The user should specify the grid when using ``-bootstrap`` or when initializing a
verification test (section :ref:`sec-verif`) or a simplified-geometry experiment (section
:ref:`sec-simp`). If one initializes PISM from a saved model state using ``-i`` then the
input file determines all grid parameters. For instance, the command

.. code-block:: none

   pismr -i foo.nc -y 100

should work fine if ``foo.nc`` is a PISM output file. Because ``-i`` input files take
precedence over options,

.. code-block:: none

   pismr -i foo.nc -Mz 201 -y 100

will give a warning that "``PISM WARNING: ignoring command-line option '-Mz'``".

.. _sec-domain-dstribution:

Parallel domain distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When running PISM in parallel with ``mpiexec -n N``, the horizontal grid is distributed
across `N` processes [#]_. PISM divides the grid into `N_x` parts in the `x` direction and
`N_y` parts in the `y` direction. By default this is done automatically, with the goal
that `N_x\times N_y = N` and `N_x` is as close to `N_y` as possible. Note that `N` should,
therefore, be a composite (not prime) number.

Users seeking to override this default can specify `N_x` and `N_y` using the :opt:`-Nx`
and :opt:`-Ny` command-line options.

Once `N_x` and `N_y` are computed, PISM computes sizes of sub-domains `M_{x,i}` so that
`\sum_{i=1}^{N_x}M_{x,i} = M_x` and `M_{x,i} - \left\lfloor M_x / N_x \right\rfloor < 1`.
To specify strip widths `M_{x,i}` and `M_{y,i}`, use command-line options :opt:`-procs_x`
and :opt:`-procs_y`. Each option takes a comma-separated list of numbers as its argument.
For example,

.. code-block:: none

   mpiexec -n 3 pisms -Mx 101 -My 101 \
                      -Nx 1 -procs_x 101 \
                      -Ny 3 -procs_y 20,61,20

splits a `101 \times 101` grid into 3 strips along the `x` axis.

To see the parallel domain decomposition from a completed run, see the :var:`rank`
variable in the output file, e.g. using ``-o_size big``. The same :var:`rank` variable is
available as a spatial diagnostic field (section :ref:`sec-saving-diagnostics`).

.. rubric:: Footnotes

.. [#] In most cases one process corresponds to one "core" of your computer.
