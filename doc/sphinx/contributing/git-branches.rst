.. include:: ../global.txt

.. default-role:: literal

.. _sec-git-branches:

Git branches
^^^^^^^^^^^^

PISM development loosely follows the Git branching model described in `A successful Git
branching model <git-branching-model_>`_ by Vincent Driessen.\ [#]_

PISM's repository contains two long-lived branches: `master` and `dev`.

The default branch `master` contains released code. This way users can clone
|pism-github-url| and get the latest PISM release.

The `dev` branch contains code that is ready to be included in the next release.

Current development is done in *topic branches* started from `dev`. Once a new feature or
improvement is finished, tested, and documented, the topic branch is merged into `dev` and
deleted. Each topic branch should contain changes related to one particular topic.\ [#]_

If one person is responsible for working on a topic branch it works well to keep it up to
date with `dev` by rebasing it on top of `dev`, effectively applying all the changes
contained in it to the current state of `dev`. If rebasing is impractical one could merge
`dev` into a topic branch to get access to some features that were not available when the
branch was started. However, merging `dev` into a topic branch "just to stay up to date"
is not a good idea since it confuses commit history.

When a released version of the code needs a fix, a "bug-fix" branch is created from the
`master` branch. When the implementation of a fix is complete, the bug-fix branch is
merged into `master` (and `master` is tagged to mark the new bug-fix release) and into
`dev` so that the fix is included in the next major release.

In |pism-github-url| branches are named using the name of the person responsible for the
branch as a prefix. For example, `ckhroulev/pnetcdf` is the name of Constantine Khroulev's
branch containing improvements of the I/O code using PnetCDF_.

.. note::

   Please commit all your changes to "topic" branches. The `master` and `dev` branches are
   managed by PISM developers at UAF.

.. rubric:: Footnotes

.. [#] This model may not be perfect but works well for a project of PISM's size.
.. [#] See `Fun with merges and purposes of branches <git-fun-with-merges_>`_ by Junio C
       Hamano for more about "topic branches."
