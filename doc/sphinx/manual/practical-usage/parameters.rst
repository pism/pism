.. include:: ../../global.txt

.. _sec-pism-defaults:

PISM's configuration parameters and how to change them
------------------------------------------------------

PISM's behavior depends on values of many flags and physical parameters (see
:ref:`sec-parameter-list` for details). Most of parameters have default values [#]_ which
are read from the configuration file |config-file| in the ``lib`` sub-directory.

It is possible to run PISM with an alternate configuration file using the :opt:`-config`
command-line option:

.. code-block:: none

   pismr -i foo.nc -y 1000 -config my_config.nc

The file ``my_config.nc`` has to contain *all* of the flags and parameters present in
|config-file|.

The list of parameters is too long to include here; please see the
:ref:`sec-parameter-list` for an automatically-generated table describing them.

Some command-line options *set* configuration parameters; some PISM executables have
special parameter defaults. To examine what parameters were used in a particular run, look
at the attributes of the ``pism_config`` variable in a PISM output file.

.. _sec-parameter-studies:

Managing parameter studies
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Keeping all PISM output files in a parameter study straight can be a challenge. If the
parameters of interest were controlled using command-line options then one can use
``ncdump -h`` and look at the ``history`` global attribute.

Alternatively, one can change parameter values by using an "overriding" configuration
file. The :opt:`-config_override` command-line option provides this alternative. A file
used with this option can have a subset of the configuration flags and parameters present
in |config-file|. Moreover, PISM adds the ``pism_config`` variable with values used in a
run to the output file, making it easy to see which parameters were used.

Here's an example. Suppose we want to compare the dynamics of an ice-sheet on Earth to the
same ice-sheet on Mars, where the only physical change was to the value of the
acceleration due to gravity. Running

.. code-block:: none

   pismr -i input.nc -y 1e5 -o earth.nc <other PISM options>

produces the "Earth" result, since PISM's defaults correspond to this planet. Next, we
create ``mars.cdl`` containing the following:

.. code-block:: none

   netcdf mars {
       variables:
       byte pism_overrides;
       pism_overrides:constants.standard_gravity = 3.728;
       pism_overrides:constants.standard_gravity_doc = "m s-2; standard gravity on Mars";
   }


Notice that the variable name is ``pism_overrides`` and not ``pism_config`` above. Now

.. code-block:: none

   ncgen -o mars_config.nc mars.cdl
   pismr -i input.nc -y 1e5 -config_override mars_config.nc -o mars.nc <other PISM options>

will create ``mars.nc``, the result of the "Mars" run. Then we can use ``ncdump`` to see
what was different about ``mars.nc``:

.. code-block:: diff

   ncdump -h earth.nc | grep pism_config: > earth_config.txt
   ncdump -h mars.nc | grep pism_config: > mars_config.txt
   diff -U 1 earth_config.txt mars_config.txt
   --- earth_config.txt	2015-05-08 12:44:43.000000000 -0800
   +++ mars_config.txt	2015-05-08 12:44:51.000000000 -0800
   @@ -734,3 +734,3 @@
                   pism_config:ssafd_relative_convergence_units = "1" ;
   -               pism_config:constants.standard_gravity_doc = "acceleration due to gravity on Earth geoid" ;
   +               pism_config:constants.standard_gravity_doc = "m s-2; standard gravity on Mars" ;
                   pism_config:constants.standard_gravity_type = "number" ;
   @@ -1057,3 +1057,3 @@
                   pism_config:ssafd_relative_convergence = 0.0001 ;
   -               pism_config:constants.standard_gravity = 9.81 ;
   +               pism_config:constants.standard_gravity = 3.728 ;
                   pism_config:start_year = 0. ;

.. _sec-saving-pism-config:

Saving PISM's configuration for post-processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to saving ``pism_config`` in the output file, PISM automatically adds this
variable to all files it writes (snap shots, time series of scalar and spatially-varying
diagnostic quantities, and backups). This may be useful for post-processing and analysis
of parameter studies as the user has easy access to all configuration options, model
choices, etc., without the need to keep run scripts around.

.. rubric:: Footnotes

.. [#] For ``pismr``, grid parameters ``Mx``, ``My``, that must be set at bootstrapping,
       are exceptions.

       .. Note: This is because we don't have a way to tell if a parameter value is a
          default *or* a conscious user choice, but we do know that command-line options
          are meant to override defaults.

          The desired behavior at bootstrapping is:

          In absence of -Mx and -My use the grid size from the input file. Values set
          with -Mx and -My override ones read from the input file. The user can set one
          (or both) -Mx and -My.
