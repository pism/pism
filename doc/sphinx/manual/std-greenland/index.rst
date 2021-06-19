.. include:: ../../global.txt

.. _sec-start:

Getting started: a Greenland ice sheet example
==============================================

This introduction is intended to be interactive and participatory, and it should work on
*your personal machine* as well as on a supercomputer. Please try the commands and view
the resulting files. Do the runs with your own values for the options. We can't hide the
fact that PISM has lots of "control knobs," but fiddling with them will help you get
going. Give it a try!

We get started with an extended example showing how to generate initial states for
prognostic model experiments on the Greenland ice sheet. Ice sheet and glacier model
studies often involve modeling present and past states using actions like the ones
demonstrated here. Our particular choices made here are motivated by the evaluation of
initialization methods in :cite:`AschwandenAdalgeirsdottirKhroulev`.

We use data assembled by the `Sea-level Response to Ice Sheet Evolution (SeaRISE)
<searise_>`_ assessment process :cite:`Bindschadler2013SeaRISE`. SeaRISE is a
community-organized assessment process that provided an upper bound on ice sheet
contributions to sea level in the next 100--200 years for the IPCC AR5 report in 2013.

This example is a hands-on first look at PISM. It is not an in-depth tutorial, and some
details of what is happening are only explained later in this Manual, which thoroughly
discusses PISM options, nontrivial modeling choices, and how to preprocess input data.

The basic runs here, mostly on coarse `20` and `10\km` grids, can be
done on a typical workstation or laptop. PISM is, however, designed to make high
resolution (e.g. `5\km` to `\sim 500\,\text{m}` grids for
whole-Greenland ice sheet modeling) possible by exploiting large-scale parallel
processing. See :cite:`AschwandenAdalgeirsdottirKhroulev`, :cite:`Golledgeetal2012`,
:cite:`Golledgeetal2013`, among other published high-resolution PISM examples.

.. toctree::

   input-data.rst

   run-1.rst

   run-1-watching.rst

   run-2.rst

   run-3.rst

   run-4.rst

   grid-sequencing.rst

   parameter-study.rst

