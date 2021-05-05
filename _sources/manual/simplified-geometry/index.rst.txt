.. include:: ../../global.txt

.. _sec-simp:

Simplified geometry experiments
===============================

There have been three stages of ice sheet model intercomparisons based on simplified
geometry experiments since the early 1990s :cite:`BuelerSpray`.

EISMINT I (European Ice Sheet Modeling INiTiative) :cite:`EISMINT96`\ [#]_ was the first of
these and involved only the isothermal shallow ice approximation (SIA). Both fixed margin
and moving margin experiments were performed in EISMINT I, and various conclusions were
drawn about the several numerical schemes used in the intercomparison. EISMINT I is
superceded, however, by verification using the full variety of known exact solutions to
the isothermal SIA :cite:`BLKCB`. The "rediscovery", since EISMINT I, of the Halfar similarity
solution with zero accumulation :cite:`Halfar83`, and verification runs using that solution,
already suffices to measure the isothermal SIA performance of PISM more precisely than
would be allowed by comparison to EISMINT I results.

EISMINT II :cite:`EISMINT00` pointed out interesting and surprising properties of the
thermocoupled SIA. References :cite:`BBL`, :cite:`Hindmarsh04`, :cite:`Hindmarsh06`,
:cite:`PayneBaldwin`, :cite:`SaitoEISMINT`, :cite:`BBssasliding` each interpret the
EISMINT II experiments and/or describe attempts to add more complete physical models to
"fix" the (perceived and real) shortfalls of ice sheet model behavior on EISMINT II
experiments. We believe that the discussion in :cite:`PayneDongelmans`,
:cite:`PayneBaldwin`, :cite:`BBL` adequately explains the "spokes" in EISMINT II
experiment F as a genuine fluid instability, while :cite:`Fowler01` and Appendix B of
:cite:`BBssasliding` adequately cautions against the continuum model that generates the
"spokes" in EISMINT II experiment H. Thus we can move on from that era of controversy. In
any case, PISM has built-in support for all of the published and unpublished EISMINT II
experiments; these are described in the next subsection.

The ISMIP (Ice Sheet Model Intercomparison Project) [#]_ round of intercomparisons covers
2008--2013 (at least). There are four components of ISMIP substantially completed, namely
HOM = Higher Order Models :cite:`ISMIPHOM`, :cite:`HOMelmer`, HEINO = Heinrich Event
INtercOmparison :cite:`GreveTakahamaCalov`, :cite:`Calovetal2009HEINOfinal`, MISMIP
(below), and MISMIP3d (also below).

PISM did not participate in ISMIP-HOM but does support most of prescribed experiments (see
:ref:`sec-ISMIP-HOM`).

PISM participated in HEINO, but this ability is unmaintained. We believe the continuum
problem described by HEINO, also used in EISMINT II experiment H (above), is not
meaningfully approximate-able because of a required discontinuous jump in the basal
velocity field. The continuum problem predicts infinite vertical velocity because of this
jump (:cite:`BBssasliding`, Appendix B). Details of the numerical schemes and their results are
irrelevant if the continuum model makes such a prediction. PISM offers the physical
continuum model described in :cite:`BBssasliding`, an SIA+SSA hybrid, as an alternative to the
continuum model used in ISMIP-HEINO and EISMINT II experiment H. Indeed the SIA+SSA hybrid
is offered as a unified shallow model for real ice sheets (section :ref:`sec-dynamics`).

A third and fourth ISMIP parts are the two parts of the Marine Ice Sheet Model
Intercomparison Project, MISMIP :cite:`MISMIP2012` and MISMIP3D :cite:`MISMIP3d2013`. These
experiments are supported in PISM, as described in subsections :ref:`sec-MISMIP` and
:ref:`sec-MISMIP3d` below.

.. toctree::
   :caption: Intercomparison projects

   eismint-2.rst

   mismip.rst

   mismip3d.rst

   ismip-hom.rst

.. rubric:: Footnotes

.. [#] See |eismint-url|

.. [#] See |ismip-url|
