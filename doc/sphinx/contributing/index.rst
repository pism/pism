.. include:: ../global.txt

.. default-role:: literal

.. _sec-contributing:

Contributing to PISM
====================

Bug reports, contributions of code, documentation, and tests are always appreciated.

You will need a GitHub_ account and some familiarity with Git_\ [#]_.

Please see :ref:`sec-bug-reports` for bug reporting guidelines.

.. note::

   By submitting code, the contributor gives irretrievable consent to the redistribution
   and modification of the contributed source code as described in the PISM's open source
   license.

Contributions are preferred via pull requests to |pism-github-url|.

#. `Fork PISM's repository. <github-help-fork_>`_
#. Create a branch that will contain your changes.
#. Implement proposed changes.

   a. Make changes to the code or documentation (or both).
   b. Test your changes.
   c. Add verification or regression tests (optional but **strongly encouraged**).
   d. Update documentation, if necessary.
   e. Update the change log `CHANGES.rst`. If your contribution contains a bug fix,
      please describe the bug and its effects.

#. `Create a pull request <github-pull-request-create_>`_ and make sure to `allow
   edits from maintainers. <github-pull-request-allow-edits_>`_

   Every time you push your code to GitHub_ CircleCI_ will

   - build it with pedantic compiler settings, treating all compiler warnings as errors
   - run `make test` in the build directory.

   Please make sure all tests pass (you will get an e-mail if there was a failure).

If you are planning a large contribution we encourage you to open an issue at
|pism-issues-url| or e-mail us at |pism-email| and interact with us frequently to ensure
that your effort is well-directed.

When working on a large contribution it is essential to stay in touch with PISM's
maintainers.

We urge you to use a version control system and give PISM maintainers access to a
repository containing your code. We realize that your group's policy may make it
impossible to put the code in a public repository. However, it is very easy to set up a
private repository instead.

.. warning::

   Working in isolation will lead to a waste of many person-hours of effort once your work
   is ready to be merged into PISM.

See sections listed below for various technical details.

.. toctree::
   :caption: Contents
   :titlesonly:

   bug-reporting.rst

   coding_guidelines.rst

   git-introduction.rst

   git-branches.rst

   development-workflow.rst

   how-to.rst

.. rubric:: Footnotes

.. [#] Please see :ref:`sec-git-introduction` for a brief introduction and `Git
       documentation`_ for more.
