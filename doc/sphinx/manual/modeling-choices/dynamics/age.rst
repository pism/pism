.. include:: ../../../global.txt

.. _sec-age:

Computing ice age
-----------------

By default, PISM does not compute the age of the ice because it does not directly impact
ice flow when using the default flow laws. It is very easy to turn on. Just set
:opt:`-age`. A 3D variable ``age`` will appear in output files. It is read at input if
``-age`` is set and otherwise it is ignored even if present in the input file. If ``-age``
is set and the variable ``age`` is absent in the input file then the initial age is set to
zero.

The age of the ice can be used in two parameterizations in the SIA stress balance model:

#. Ice grain size parameterization based on data from :cite:`DeLaChapelleEtAl98` and
   :cite:`LipenkovEtAl89` (Vostok core data). In PISM, only the Goldsby-Kohlstedt flow law
   (see :ref:`sec-rheology`) uses the grain size.

#. The flow enhancement factor can be coupled to the age of the ice as in
   :cite:`Greve97Greenland`: during Eemian and Holocene `e` is set to

   - :config:`stress_balance.sia.enhancement_factor_interglacial` during Eemian and Holocene,
   - :config:`stress_balance.sia.enhancement_factor` otherwise.

See :config:`time.eemian_start`, :config:`time.eemian_end`, and
:config:`time.holocene_start`.
