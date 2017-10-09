.. include:: ../global.txt

.. _sec-users-manual:

PISM User's Manual
==================

Welcome!  All information about PISM is online at the home page

    |pism-docs|

Please see the extensive lists of PISM publications and applications at that page.

This User's Manual gives examples of how to run PISM using publicly-available data for:
the whole Greenland ice sheet, the Jakobshavn outlet glacier in Greenland, the Ross ice
shelf in Antarctica, and a number of simplified geometry tests. It documents all the
PISM's command-line options and configuration parameters. It summarizes the continuum
models used by PISM, and it illustrates how PISM's numerical approximations are verified.

See the :ref:`Installation Manual <sec-installation>` for how to download the PISM source
code and install it, along with needed libraries. The :ref:`Climate Forcing Manual
<sec-climate-forcing>` extends this Manual to cover additional couplings to atmosphere and
ocean models and data.

Users who want to understand more deeply how PISM is designed, or who want to extend it,
will need to go beyond what is described here. See the `Source Code Browser
<pism-browser_>`_, which is online for the latest stable version. It can be generated from
source code as described in :ref:`sec-install-documentation`. It gives a complete view of the
class/object structure of the PISM source code.

.. warning::

   PISM is an ongoing research project. Ice sheet modeling requires many choices. Please
   don't trust the results of PISM or any other ice sheet model without a fair amount of
   exploration. Also, please don't expect all your questions to be answered here. `Write to
   us <pism-email_>`_ with questions.

.. toctree::
   :caption: Contents

   std-greenland/index.rst

   highlevelview/index.rst

   initialization/index.rst

   modeling-choices/index.rst

   practical-usage/index.rst

   simplified-geometry/index.rst

   verification/index.rst

   validation/index.rst

   jakobshavn/index.rst

   parameters/index.rst

   diagnostics/index.rst
