PISM release checklist
======================

#. [ ] Create a new ``stableXX`` branch, if needed.
#. [ ] Merge, if needed.
#. [ ] Set ``Pism_BRANCH`` in ``CMakeLists.txt``
#. [ ] Replace ``0.x`` with ``0.x+1`` in docs.
#. [ ] Update ``CHANGES.rst``.
#. [ ] Tag.

   ::

       git tag -a v0.X -m "The v0.X release. See CHANGES.rst for the list of changes since v0.X-1."

#. [ ] Push.

   ::

       git push -u origin HEAD

#. [ ] Push tags.

   ::

       git push --tags

#. [ ] Re-build docs.

   ::

       make manual_html manual_pdf browser.tgz pismpython_docs

#. [ ] Upload these docs.
#. [ ] Switch the default branch on https://github.com/pism/pism
#. [ ] Write a news item for ``pism-docs.org``.
#. [ ] Replace ``0.x`` with ``0.x+1`` on ``pism-docs.org``.
#. [ ] Send an e-mail to CRYOLIST.
#. [ ] Tell more people, if desired.
