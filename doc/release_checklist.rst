PISM release checklist
======================

#. Run ``make`` in the ``doc/sphinx`` directory to update lists of diagnostics and
   configuration parameters.
#. Run ``make`` in the ``doc`` directory to update funding sources.
#. Create a new branch and merge (if needed).
#. Set ``Pism_BRANCH`` in ``CMakeLists.txt`` to "stable".
#. Update ``version`` and ``release`` in ``doc/sphinx/conf.py``.
#. Update ``CHANGES.rst``.
#. Tag.

   ::

      git tag -a v0.X -m "The v0.X release. See CHANGES.rst for the list of changes since v0.X-1."

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
#. Update the default branch on https://github.com/pism/pism
#. Write a news item for ``pism-docs.org``.
#. Update the current PISM version on ``pism-docs.org``.
#. Send an e-mail to CRYOLIST.
#. Tell more people, if desired.
