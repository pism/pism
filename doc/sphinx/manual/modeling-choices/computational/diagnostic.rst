.. include:: ../../../global.rst

.. _sec-diagnostic-computations:

Diagnostic computations
-----------------------

A "diagnostic" computation can be defined as one where the internal state does not evolve.
The internal state of PISM is the set of variables read by "``-i``". You can ask PISM to
do a diagnostic computation by setting the run duration to a small number such as
`0.001` years (about `9` hours). The duration to use depends on the modeling
setup, but should be smaller than the maximum time-step allowed by PISM's stability
criteria. Such short runs can also be used to look at additional fields corresponding to
the current model state.

As an example, consider these two runs:

.. code-block:: none

   pisms -y 6000 -o foo.nc
   pismr -i foo.nc -y 0.001 -o bar.nc -o_size big

The result of the second (short) run is a NetCDF file ``bar.nc`` which contains the full
three-dimensional velocity field in the scalar NetCDF variables ``uvel``, ``vvel``, and
``wvel``, as well as many other variables. The file ``foo.nc`` does not contain many of
these fields because it was written with the default output size of ``medium``. The "``-y
0.001``" run has diagnostically "filled-in" all the fields which PISM can model at a time
step, but the run duration was chosen so as to avoid significant model state evolution
during the run.

This diagnostic mode is often associated to the modeling of ice shelves and ice streams.
section :ref:`sec-ross` describes using a short "diagnostic" run to model the Ross ice
shelf :cite:`MacAyealetal`. Verification tests I and J, section :ref:`sec-verif`, are
diagnostic calculations using the SSA.

The NetCDF model state saved by PISM at the end of an *evolution* run (i.e. with "``-y
Y``" for `Y>0`) does not, under the default ``-o_size medium`` output size, contain
the three-dimensional velocity field. Instead, it contains just a few more variables than
those which are needed to restart the run with ``-i``. One can force PISM to save all the
supported diagnostic quantities at the end of a time-stepping run using the option
``-o_size big``. Or one can go back and do a "``-y small_number``" diagnostic run using
``-o_size big``.
