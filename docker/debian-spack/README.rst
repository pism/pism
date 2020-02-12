This directory contains a `Dockerfile` that can be used to test PISM's installation with
Spack.

Run

.. code-block:: bash

   docker build -t some-tag .

to try installing PISM on Debian 10.

The file `packages.yaml` tells Spack to use Debian's binaries for some of PISM's
dependencies.
