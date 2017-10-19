.. include:: ../../global.txt

.. _sec-validation:

Validation case studies
=======================


"Validation" describes the comparison of numerical model output with physical observations
in cases where the observations are sufficiently-complete and of sufficient quality so
that the performance of the numerical model can be assessed :cite:`Roache`, :cite:`Wesseling`.
Roughly speaking, validation can happen when the observations or data are better than the
model, so the comparison measures the quality of the numerical model and not merely errors
in, or incompleteness of, the data. Because of the difficulty of measuring boundary
conditions for real ice flows, this situation is not automatic in glaciology, or even
common.\ [#]_ Nonetheless we try two cases, first where PISM is applied on
millimeter scale to model a laboratory experiment, and second for a large-scale ice flow
in which all uncertainties of bedrock topography, basal sliding, and subglacial hydrology
are removed, namely a present-day ice shelf.

.. toctree::

   labgum.rst

   ross.rst

.. rubric:: Footnotes

.. [#] Which explains the rise of "simplified geometry intercomparisons"; see section
       :ref:`sec-simp`.
