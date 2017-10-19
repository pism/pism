.. include:: ../../../global.txt

.. _sec-bedsmooth:

Parameterization of bed roughness in the SIA
--------------------------------------------

Schoof :cite:`Schoofbasaltopg2003` describes how to alter the SIA stress balance to model
ice flow over bumpy bedrock topgraphy. One computes the amount by which bumpy topography
lowers the SIA diffusivity. An internal quantity used in this method is a smoothed version
of the bedrock topography. As a practical matter for PISM, this theory improves the SIA's
ability to handle bed roughness because it parameterizes the effects of "higher-order"
stresses which act on the ice as it flows over bumps. For additional technical description
of PISM's implementation, see :ref:`sec-bed-roughness`.

There is only one associated option: :opt:`-bed_smoother_range` gives the half-width of
the square smoothing domain in meters. If zero is given, ``-bed_smoother_range 0`` then
the mechanism is turned off. The mechanism is on by default using executable ``pismr``,
with the half-width set to 5 km (``-bed_smoother_range 5.0e3``), giving Schoof's
recommended smoothing size of 10 km :cite:`Schoofbasaltopg2003`.

This mechanism is turned off by default in executables ``pisms`` and ``pismv``.

Under the default setting ``-o_size medium``, PISM writes fields :var:`topgsmooth` and
:var:`schoofs_theta` from this mechanism. The thickness relative to the smoothed bedrock
elevation, namely :var:`topgsmooth`, is the difference between the unsmoothed surface
elevation and the smoothed bedrock elevation. It is *only used internally by this
mechanism*, to compute a modified value of the diffusivity; the rest of PISM does not use
this or any other smoothed bed. The field :var:`schoofs_theta` is a number `\theta`
between `0` and `1`, with values significantly below zero indicating a reduction in
diffusivity, essentially a drag coefficient, from bumpy bed.
