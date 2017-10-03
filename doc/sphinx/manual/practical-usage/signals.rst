.. include:: ../../global.txt

.. |pid| replace:: *PID*\s

.. _sec-signal:

Signals, to control a running PISM model
----------------------------------------

Ice sheet model runs sometimes take a long time, so the state of a run may need checking.
Sometimes the run needs to be stopped, but with the possibility of restarting. PISM
implements these behaviors using "signals" from the POSIX standard, included in Linux and
most flavors of Unix. :numref:`tab-signals` summarizes how PISM responds to signals.
A convenient form of ``kill``, for Linux users, is ``pkill`` which will find processes by
executable name. Thus "``pkill -USR1 pismr``" might be used to send all PISM processes the
same signal, avoiding an explicit list of |pid|.

.. list-table:: Signalling running PISM processes.  "|pid|" stands for list of all identifiers of the PISM processes.
   :name: tab-signals
   :header-rows: 1
   :widths: 1,1,2

   * - Command
     - Signal
     - PISM behavior
   * - ``kill -KILL`` |pid|
     - ``SIGKILL``
     - Terminate with extreme prejudice. PISM cannot catch it and no state is saved.
   * - ``kill -TERM`` |pid|
     - ``SIGTERM``
     - End process(es), but save the last model state in the output file, using ``-o``
       name or default name as normal. Note that the ``history`` string in the output file
       will contain an "``EARLY EXIT caused by signal SIGTERM``" indication.
   * - ``kill -USR1`` |pid|
     - ``SIGUSR1``
     - Process(es) will continue after saving the model state at the end of the current
       time step, using a file name including the current model year. Time-stepping is not
       altered. Also flushes output buffers of scalar time-series.
   * - ``kill -USR2`` |pid|
     - ``SIGUSR2``
     - Just flush time-series output buffers.
   
Here is an example. Suppose we start a long verification run in the background, with
standard out redirected into a file:

.. code-block:: none

   pismv -test G -Mz 101 -y 1e6 -o testGmillion.nc >> log.txt &

This run gets a Unix process id, which we assume is "8920". (Get it using ``ps`` or
``pgrep``.) If we want to observe the run without stopping it we send the ``USR1`` signal:


.. code-block:: none

   kill -USR1 8920

(With ``pkill`` one can usually type "``pkill -usr1 pismv``".) Suppose it happens that we
caught the run at year 31871.5. Then, for example, a NetCDF file ``pismv-31871.495.nc`` is
produced. Note also that in the standard out log file ``log.txt`` the line

.. code-block:: none

   caught signal SIGUSR1:  Writing intermediate file ... and flushing time series.

appears around that time step. Suppose, on the other hand, that the run needs to be
stopped. Then a graceful way is

.. code-block:: none

   kill -TERM 8920

because the model state is saved and can be inspected.
