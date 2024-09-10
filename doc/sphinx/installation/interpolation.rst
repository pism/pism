.. include:: ../global.txt

.. default-role:: literal

.. _sec-install-yac:

Enabling flexible interpolation using YAC and PROJ
--------------------------------------------------

PISM's projection-aware interpolation code uses PROJ_\ [#f1]_ to compute *longitude, latitude*
coordinates of cell centers and cell corners and YAC_ for the interpolation itself.

YAC depends on YAXT_. Use the following commands to install it:

.. literalinclude:: code/yaxt.sh
   :language: bash
   :caption: Building YAXT
   :start-after: manual-begin
   :end-before: manual-end

Once this is done, build and install YAC:

.. literalinclude:: code/yac.sh
   :language: bash
   :caption: Building YAC
   :start-after: manual-begin
   :end-before: manual-end

To build PISM with these libraries, set the CMake variable `Pism_USE_YAC_INTERPOLATION`
and include YAC's installation prefix in `CMAKE_PREFIX_PATH`:

.. code-block:: bash

   cmake -DPism_USE_YAC_INTERPOLATION=YES \
         -DCMAKE_PREFIX_PATH="$HOME/local/yac;..." \
         ...

See section :ref:`sec-interpolation` to learn how to use flexible interpolation in PISM.

.. rubric:: Footnotes

.. [#f1] PROJ is widely used and you should not have to build it from sources.
