.. include:: ../../global.txt

.. _sec-snapshots:

Snapshots of the model state
----------------------------

Sometimes you want to check the model state every `1000` years, for example. One possible
solution is to run PISM for a thousand years, have it save all the fields at the end of
the run, then restart and run for another thousand, and etc. This forces the adaptive
time-stepping mechanism to stop *exactly* at multiples of `1000` years, which may be
desirable in some cases.

If saving exactly at specified times is not critical, then use the :opt:`-save_file` and
:opt:`-save_times` options. For example,

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -save_file snapshots.nc \
         -save_times 1000:1000:10000

starts a PISM evolution run, initializing from ``foo.nc``, running for 10000 years and
saving snapshots to ``snapshots.nc`` at the first time-step after each of the years `1000`,
`2000`, ..., `10000`.

We use a MATLAB-style range specification, `a:\dt:b`, where `a,\dt,b` are in years. The
time-stepping scheme is not affected, but as a consequence we do not guarantee producing
the exact number of snapshots requested if the requested save times have spacing
comparable to the model time-steps. This is not a problem in the typical case in which
snapshot spacing is much greater than the length of a typical time step.

It is also possible to save snapshots at intervals that are not equally-spaced by giving
the :opt:`-save_times` option a comma-separated list. For example,

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc \
         -save_file snapshots.nc \
         -save_times 1000,1500,2000,5000

will save snapshots on the first time-step after years `1000`, `1500`, `2000` and `5000`.
The comma-separated list given to the :opt:`-save_times` option can be at most `200`
numbers long.

If ``snapshots.nc`` was created by the command above, running

.. code-block:: none

   pismr -i snapshots.nc -y 1000 -o output_2.nc

will initialize using the last record in the file, at about :math:`5000` years. By
contrast, to restart from :math:`1500` years (for example) it is necessary to extract the
corresponding record using ``ncks``

.. code-block:: none

   ncks -d t,1500years snapshots.nc foo.nc

and then restart from ``foo.nc``. Note that ``-d t,N`` means "extract the :math:`N`-th
record" (counting from zero). So, this command is equivalent to

.. code-block:: none

   ncks -d t,1 snapshots.nc foo.nc

Also note that the second snapshot will probably be *around* :math:`1500` years and
``ncks`` handles this correctly: it takes the record closest to :math:`1500` years.

By default re-startable snapshots contain only the variables needed for restarting PISM.
Use the command-line option :opt:`-save_size` to change what is saved.

Another possible use of snapshots is for restarting runs on a batch system which kills
jobs which go over their allotted time. Running PISM with options ``-y 1500``
``-save_times 1000:100:1400`` would mean that if the job is killed before completing the
whole 1500 year run, we can restart from near the last multiple of `100` years.
Restarting with option :opt:`-ye` would finish the run on the desired year.

When running PISM on such a batch system it can also be useful to save re-startable
snapshots at equal wall-clock time (as opposed to model time) intervals by adding the
":opt:`-backup_interval` (hours)" option.

.. caution::

   If the wall-clock limit is equal to `N` times backup interval for a whole number
   `N` PISM will likely get killed while writing the last backup.

It is also possible to save snapshots to separate files using the :opt:`-save_split` option.
For example, the run above can be changed to

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -save_file snapshots \
         -save_times 1000,1500,2000,5000 -save_split

for this purpose. This will produce files called ``snapshots-year.nc``. This option is
generally faster if many snapshots are needed, apparently because of the time necessary to
reopen a large file at each snapshot when :opt:`-save_split` is not used. Note that tools
like NCO and ``ncview`` usually behave as desired with wildcards like
"``snapshots-*.nc``".

.. _sec-snapshot-parameters:

Parameters
==========

.. pism-parameters::
   :prefix: output.snapshot.
