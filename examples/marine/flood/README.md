### A simplified setup showing that `-kill_icebergs` should always be "on".

We start with three blobs of ice on a flat bed. As the run goes on,
the sea level rises, and blobs start floating, starting with the
smallest one.

Once the whole blob is afloat, the the SSA solver fails because of a
zero pivot.

Note that without the calving front boundary condition (`make
no_cfbc`), the strength extension "ties" floating blobs to the
biggest, blob, which remains grounded.
