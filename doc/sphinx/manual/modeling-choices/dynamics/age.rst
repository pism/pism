.. include:: ../../../global.rst

.. _sec-age:

Computing ice age
-----------------

By default, PISM does not compute the age of the ice because it does not directly impact
ice flow when using the default flow laws. It is very easy to turn on. Just set
:opt:`-age`. A 3D variable ``age`` will appear in output files. It is read at input if
``-age`` is set and otherwise it is ignored even if present in the input file. If ``-age``
is set and the variable ``age`` is absent in the input file then the initial age is set to
zero.

..
   describe coupling of ice age and the SIA flow
