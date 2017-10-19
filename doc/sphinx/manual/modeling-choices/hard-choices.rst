.. include:: ../../global.txt

.. _sec-hard-choices:

Dealing with more difficult modeling choices
--------------------------------------------

Most uses of an ice sheet model depend on careful modeling choices in situations where
there are considerable uncertainties *and* the model results depend strongly on those
choices. There may be, at the present state of knowledge, *no clear default values* that
PISM can provide. Furthermore, the available PISM options and sub-models are known to
*not* be sufficient for all users. Thus there are modelling situations for which we know
the user may have to do a great deal more hard work than just choose among PISM runtime
options.

Here are example cases where users have worked hard:

- User made use of available data in order to choose parameters for existing PISM models.
  These parameters then override PISM defaults.

  .. admonition:: Example
     :class: note

     Use regional atmosphere model output to identify PDD parameters suitable for modeling
     surface mass balance on a particular ice sheet. Then supply these parameters to PISM
     by a ``-config_override`` file.

     .. our UAF current situation with Greenland

- User wrote code, including code which modified current PISM internals, either to add
  additional processes or to "correct" PISM default process models.

  .. admonition:: Example
     :class: note

     Add a new sub-ice-shelf melt model by modifying C++ code in the ``src/coupler/``
     directory.

     .. PIK ocean models

- User simplified the model in use, instead of the default which was more elaborate.

  .. admonition:: Example
     :class: note

     Instead of using the PISM default mechanism connecting basal melt rate and basal
     strength, bypass this mechanism by generating a map of yield stress ``tauc`` directly
     and supplying it as input.

     .. Nick's -yield_stress constant choice
