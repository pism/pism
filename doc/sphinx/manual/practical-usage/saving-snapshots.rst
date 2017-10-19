.. include:: ../../global.txt

.. _sec-snapshots:

Saving re-startable snapshots of the model state
------------------------------------------------

Sometimes you want to check the model state every 1000 years, for example. One possible
solution is to run PISM for a thousand years, have it save all the fields at the end of
the run, then restart and run for another thousand, and etc. This forces the adaptive
time-stepping mechanism to stop *exactly* at multiples of 1000 years, which may be
desirable in some cases.

If saving exactly at specified times is not critical, then use the ``-save_file`` and
``-save_times`` options. For example,

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -save_file snapshots.nc \
         -save_times 1000:1000:10000

starts a PISM evolution run, initializing from ``foo.nc``, running for 10000 years and
saving snapshots to ``snapshots.nc`` at the first time-step after each of the years 1000,
2000, ..., 10000.

We use a MATLAB-style range specification, :math:`a:\Delta t:b`, where :math:`a,\Delta
t,b` are in years. The time-stepping scheme is not affected, but as a consequence we do
not guarantee producing the exact number of snapshots requested if the requested save
times have spacing comparable to the model time-steps. This is not a problem in the
typical case in which snapshot spacing is much greater than the length of a typical time
step.

It is also possible to save snapshots at intervals that are not equally-spaced by giving
the ``-save_times`` option a comma-separated list. For example,

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -save_file snapshots.nc \
         -save_times 1000,1500,2000,5000

will save snapshots on the first time-step after years 1000, 1500, 2000 and 5000. The
comma-separated list given to the ``-save_times`` option can be at most 200 numbers long.

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
Use the command-line option ``-save_size`` to change what is saved.

Another possible use of snapshots is for restarting runs on a batch system which kills
jobs which go over their allotted time. Running PISM with options ``-y 1500``
``-save_times 1000:100:1400`` would mean that if the job is killed before completing the
whole 1500 year run, we can restart from near the last multiple of :math:`100` years.
Restarting with option ``-ye`` would finish the run on the desired year.

When running PISM on such a batch system it can also be useful to save re-startable
snapshots at equal wall-clock time (as opposed to model time) intervals by adding the
":opt:`-backup_interval` (hours)" option.

.. caution::

   If the wall-clock limit is equal to :math:`N` times backup interval for a whole number
   :math:`N` PISM will likely get killed while writing the last backup.

It is also possible to save snapshots to separate files using the ``-save_split`` option.
For example, the run above can be changed to

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -save_file snapshots \
         -save_times 1000,1500,2000,5000 -save_split

for this purpose. This will produce files called ``snapshots-year.nc``. This option is
generally faster if many snapshots are needed, apparently because of the time necessary to
reopen a large file at each snapshot when ``-save_split`` is not used. Note that tools
like NCO and ``ncview`` usually behave as desired with wildcards like
"``snapshots-*.nc``".

:numref:`tab-snapshot-opts` lists the options related to saving snapshots of the
model state.

.. list-table:: Command-line options controlling saving snapshots of the model state.
   :name: tab-snapshot-opts
   :header-rows: 1
   :widths: 1,1

   * - Option
     - Description
   * - :opt:`-save_file`
     - Specifies the file to save to.
   * - :opt:`-save_times`
     - Specifies times at which to save snapshots, by either a MATLAB-style range
       :math:`a:\Delta t:b` or a comma-separated list.
   * - :opt:`-save_split`
     - Separate the snapshot output into files named ``snapshots-year.nc``. Faster if you
       are saving more than a dozen or so snapshots.
   * - :opt:`-save_size` ``[none,small,medium,big,big_2d]``
     - Similar to ``o_size``, changes the "size" of the file (or files) written; the
       default is "small"
