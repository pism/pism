.. include:: ../../global.txt

.. _sec-diagnostic-viewers:

Run-time diagnostic viewers
---------------------------

Basic graphical views of the changing state of a PISM ice model are available at the
command line by using options listed in :numref:`tab-diag-viewers`. All the quantities
listed in :ref:`sec-extra_vars` are available. Additionally, a couple of diagnostic
quantities are *only* available as run-time viewers; these are shown in table
:numref:`tab-special-diag-viewers`.

.. list-table:: Options controlling run-time diagnostic viewers
   :name: tab-diag-viewers
   :header-rows: 1
   :widths: 1,3

   * - Option
     - Description
   * - :opt:`-view`
     - Turns on map-plane views of one or several variables, see :ref:`sec-extra_vars`.
   * - :opt:`-view_size` (number)
     - desired viewer size, in pixels
   * - :opt:`-display`
     - The option ``-display :0`` seems to frequently be needed to let PETSc use Xwindows
       when running multiple processes. It must be given as a *final* option, after all
       the others.

The option ``-view`` shows map-plane views of 2D fields and surface and basal views of 3D
fields (see :ref:`sec-extra_vars`); for example:

.. code-block:: none

   pismr -i input.nc -y 1000 -o output.nc -view thk,tempsurf

shows ice thickness and ice temperature at the surface.

.. list-table:: Special run-time-only diagnostic viewers
   :name: tab-special-diag-viewers
   :header-rows: 1
   :widths: 3,5

   * - Option
     - Description
   * - :opt:`-ssa_view_nuh`
     - log base ten of ``nuH``, only available if the finite-difference SSA solver is
       active.
   * - :opt:`-ssa_nuh_viewer_size` (number)
     - Adjust the viewer size.
   * - :opt:`-ssafd_ksp_monitor_draw`
     - Iteration monitor for the Krylov subspace routines (KSP) in PETSc. Residual norm
       versus iteration number.
