.. include:: ../global.txt

.. _sec-release-checklist:

Release checklist
=================

#. Run ``make manual_linkcheck`` and fix any broken links in the manual.
#. Run ``make`` in the ``doc/sphinx`` directory to update lists of diagnostics and
   configuration parameters.
#. Run ``make`` in the ``doc`` directory to update funding sources.
#. Create a "pre-release" branch starting from the "``dev``" branch and remove code that
   should not be a part of the release.
#. Set ``Pism_BRANCH`` in ``CMakeLists.txt`` to "``stable``".
#. Update ``version``, ``release``, and ``copyright`` in ``doc/sphinx/conf.py``.
#. Update ``CHANGES.rst``.
#. Tag.

   ::

      git tag -a v1.X -m "The v1.X release. See CHANGES.rst for the list of changes since v1.X-1."

#. Push.

   ::

      git push -u origin HEAD

#. Push tags.

   ::

      git push --tags

#. Re-build docs.

   ::

      make manual_html manual_pdf browser.tgz pismpython_docs

#. Upload these docs.
#. Write a news item for ``pism-docs.org``.
#. Update the current PISM version on ``pism-docs.org``.
#. Send an e-mail to CRYOLIST.
#. Tell more people, if desired.
#. Create a new "release" on https://github.com/pism/pism/releases
