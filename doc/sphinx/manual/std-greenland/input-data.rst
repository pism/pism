.. include:: ../../global.txt

Input data
----------

The NetCDF data used to initialize SeaRISE runs is `freely-available online
<searise-greenland_>`_.

To download the specific file we want, namely ``Greenland_5km_v1.1.nc``, and preprocess it
for PISM, do:

.. code-block:: none

   cd examples/std-greenland
   ./preprocess.sh

The script ``preprocess.sh`` requires ``wget`` and also the `NetCDF Operators <NCO_>`_. It
downloads the version 1.1 of the SeaRISE "master" present-day data set, which contains ice
thickness and bedrock topography from BEDMAP :cite:`BamberLayberryGogenini`, and modeled
precipitation and surface mass balance rates from RACMO :cite:`Ettemaetal2009`, among
other fields.

In particular, it creates three new NetCDF files which can be read by PISM. The
spatially-varying fields, with adjusted metadata, go in ``pism_Greenland_5km_v1.1.nc``.
The other two new files contain famous time-dependent paleo-climate records from ice and
seabed cores: ``pism_dT.nc`` has the GRIP temperature record :cite:`JohnsenetalGRIP` and
``pism_dSL.nc`` has the SPECMAP sea level record :cite:`Imbrieetal1984`.

Any of these NetCDF files can be viewed with ``ncview`` or other NetCDF visualization
tools; see :numref:`tab-NetCDFview`. An application of IDV to the master data set produced
:numref:`fig-sr-input`, for example. Use ``ncdump -h`` to see the metadata and history of
the files.

.. figure:: figures/sr-greenland-input.png
   :name: fig-sr-input

   The input file contains present-day ice thickness (left; m), bedrock elevation (center;
   m), and present-day precipitation (right; `m / year` ice equivalent) for
   SeaRISE-Greenland. These are fields :var:`thk`, :var:`topg`, and :var:`precipitation`,
   respectively, in ``pism_Greenland_5km_v1.1.nc``.
