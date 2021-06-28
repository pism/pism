.. include:: ../../global.txt

.. _sec-EISMINTII:

EISMINT II
----------

There are seven experiments described in the published EISMINT II writeup :cite:`EISMINT00`.
They are named A, B, C, D, F, G, and H. They have these common features:

- runs are for `2\times 10^5` years, with no prescribed time step;
- a `61\times 61` horizontal grid on a square domain (`1500` km side length) is prescribed;
- surface inputs (temperature and mass balance) have angular symmetry around the grid center;
- the bed is flat and does not move (no isostasy);
- the temperature in the bedrock is not modeled;
- only the cold (not polythermal) thermomechanically-coupled SIA is used :cite:`EISMINT00`; and
- basal melt rates do not affect the evolution of the ice sheet.

The experiments differ from each other in their various combinations of surface
temperature and mass balance parameterizations. Experiments H and G involve basal sliding,
under the physically-dubious SIA sliding rubric (:cite:`BBssasliding`, Appendix B), while the
others don't. Four experiments start with zero ice (A,F,G,H), while the other experiments
(B,C,D) start from the final state of experiment A.

In addition to the seven experiments published in :cite:`EISMINT00`, there were an additional
five experiments described in the EISMINT II intercomparison description :cite:`EISIIdescribe`,
labeled E, I, J, K, and L. These experiments share most features listed above, but with
the following differences. Experiment E is the same as experiment A except that the peak
of the accumulation, and also the low point of the surface temperature, are shifted by 100
km in both `x` and `y` directions; also experiment E starts with the final state of
experiment A. Experiments I and J are similar to experiment A but with non-flat "trough"
bed topography. Experiments K and L are similar to experiment C but with non-flat "mound"
bed topography.

See :numref:`tab-eisII` for how to run supported EISMINT II experiments in PISM.
Experiments E -- L are only documented in :cite:`EISIIdescribe`.

.. note::

   Experiments G and H are not supported.

.. list-table:: Running the EISMINT II experiments in PISM. Use ``-skip -skip_max 5``, on
                the `61\times 61` default grid, for significant speedup.
   :name: tab-eisII
   :header-rows: 1
   :widths: 5,2

   * - Command: "``pismr +``"
     - Relation to experiment A

   * - ``-eisII A -Mx 61 -My 61 -Mz 61 -Lz 5000 -y 2e5 -o eisIIA.nc``
     -
   * - ``-eisII B -i eisIIA.nc -y 2e5 -o eisIIB.nc``
     - warmer
   * - ``-eisII C -i eisIIA.nc -y 2e5 -o eisIIC.nc``
     - less snow (lower accumulation)
   * - ``-eisII D -i eisIIA.nc -y 2e5 -o eisIID.nc``
     - smaller area of accumulation
   * - ``-eisII F -Mx 61 -My 61 -Mz 81 -Lz 6000 -y 2e5 -o eisIIF.nc``
     - colder; famous spokes :cite:`BBL`
   * - ``-eisII E -i eisIIA.nc -y 2e5 -o eisIIE.nc``
     - shifted climate maps
   * - ``-eisII I -Mx 61 -My 61 -Mz 61 -Lz 5000 -y 2e5 -o eisIII.nc``
     - trough topography
   * - ``-eisII J -i eisIII.nc -y 2e5 -o eisIIJ.nc``
     - trough topography and less snow
   * - ``-eisII K -Mx 61 -My 61 -Mz 61 -Lz 5000 -y 2e5 -o eisIIK.nc``
     - mound topography
   * - ``-eisII L -i eisIIK.nc -y 2e5 -o eisIIL.nc``
     - mound topography and less snow

The vertical grid is not specified in EISMINT II, but a good simulation of the
thermomechanically-coupled conditions near the base of the ice requires relatively-fine
resolution there. We suggest using the default unequally-spaced grid. With 61 levels it
gives a grid spacing of `\sim 20 m` in the ice layer closest to the bed, but more vertical
levels are generally better. Alternatively these experiments can be done with an
equally-spaced grid; in this case we suggest using enough vertical levels to give 20 m
spacing, for example. When there is sliding, even more vertical resolution is recommended
(see :numref:`tab-eisII`). Also, the vertical extent must be sufficient so that when
the ice thickness grows large, especially before thermo-softening brings it back down, the
vertical grid is tall enough to include all the ice. :numref:`tab-eisII` therefore
includes suggested settings of ``-Lz``; experiment F is different because ice thickness
increases with colder temperatures.

These SIA-only simulations parallelize well. Very roughly, for the standard `61\times 61`
horizontal grid, wall-clock-time speedups will occur up to about 30 processors. Runs on
finer (horizontal) grids will benefit from even more processors. Also, the "skip"
mechanism which avoids updating the temperature at each time step is effective, so options
like ``-skip -skip_max 5`` are recommended.

The EISMINT II experiments can be run with various modifications of the default settings.
For instance, a twice-finer grid in the horizontal is "``-Mx 121 -My 121``".
:numref:`tab-eisIIoptions` lists some optional settings which are particular to the
EISMINT II experiments.

.. list-table:: Changing the default settings for EISMINT II
   :name: tab-eisIIoptions
   :header-rows: 1
   :widths: 1,2,1,2

   * - Option
     - Default values [experiments]
     - Units
     - Meaning

   * - :opt:`-eisII`
     - A
     -
     - Choose single character name of EISMINT II :cite:`EISMINT00` simplified geometry
       experiment. See :numref:`tab-eisII`.

   * - :opt:`-Mmax`
     - 0.5 [ABDEFIK],

       0.25 [CJL]
     - `m / year`
     - max value of accumulation rate

   * - :opt:`-Rel`
     - 450 [ABEFIK],

       425 [CDJL]
     - km
     - radial distance to equilibrium line

   * - :opt:`-Sb`
     - `10^{-2}` [*all*]
     - `(m/year)/km`
     - radial gradient of accumulation rate

   * - :opt:`-ST`
     - `1.67 \times 10^{-2}` [*all*]
     - K/km
     - radial gradient of surface temperature

   * - :opt:`-Tmin`
     - 238.15 [ACDEIJKL],
       
       243.15 [B],
       
       223.15 [F]
     - K
     - max of surface temperature

   * - :opt:`-bmr_in_cont`
     -
     -
     - Include the basal melt rate in the mass continuity computation; overrides EISMINT
       II default.

See subdirectory ``examples/eismintII/`` for a simple helper script ``runexp.sh``.
