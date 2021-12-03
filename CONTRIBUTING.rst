Contributing to PISM
====================

There are many ways you can contribute to PISM:

- Fix typos, inaccuracies, and omissions in the manual.
- Improve documentation of existing features.
- Provide additional examples.
- Add new tests for existing code.
- Report issues with the code or documentation.
- Fix bugs in PISM.
- Implement new features.

Please see `Contributing to PISM <pism-contributing_>`_ in PISM's manual for some guidelines.

In summary: documentation and code contributions are preferred via pull requests to
|pism-github-url|.

#. `Fork PISM's repository. <github-help-fork_>`_
#. Create a branch that will contain your changes.
#. Implement proposed changes.

   a. Make changes to the code or documentation (or both).
   b. Test your changes.
   c. Add verification or regression tests (optional but **strongly encouraged**).
   d. Update documentation, if necessary.
   e. Update the change log ``CHANGES.rst``. If your contribution contains a bug fix,
      please describe the bug and its effects.

#. `Create a pull request <github-pull-request-create_>`_ and make sure to `allow
   edits from maintainers. <github-pull-request-allow-edits_>`_

If you are planning a large contribution we encourage you to open an issue at
|pism-issues-url| or e-mail us at |pism-email| and interact with us frequently to ensure
that your effort is well-directed.

.. note::

   By submitting code, the contributor gives irretrievable consent to the redistribution
   and modification of the contributed source code as described in the PISM's open source
   license.

.. URLs

.. _github-help-fork: https://help.github.com/en/articles/fork-a-repo
.. _github-pull-request-create: https://help.github.com/en/articles/creating-a-pull-request
.. _github-pull-request-allow-edits: https://help.github.com/en/articles/allowing-changes-to-a-pull-request-branch-created-from-a-fork
.. _pism-contributing: http://www.pism.io/docs/contributing/

.. |pism-github-url| replace:: https://github.com/pism/pism
.. |pism-issues-url| replace:: https://github.com/pism/pism/issues
.. |pism-email| replace:: uaf-pism@alaska.educ

..
   Local Variables:
   fill-column: 90
   End:
