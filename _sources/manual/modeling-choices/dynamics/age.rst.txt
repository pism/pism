.. include:: ../../../global.txt

.. _sec-age:

Computing ice age
-----------------

By default, PISM does not compute the age of the ice because it does not directly impact
ice flow when using the default flow laws. Set :config:`age.enabled` to turn it on. A 3D
variable ``age`` will appear in output files. It is read at input if :config:`age.enabled`
is set and otherwise it is ignored even if present in the input file. If
:config:`age.enabled` is set and the variable ``age`` is absent in the input file then the
initial age is set to :config:`age.initial_value`.
