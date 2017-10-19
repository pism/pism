.. include:: ../../global.txt

.. _sec-code-modifications:

Managing source code modifications
----------------------------------

"Practical usage" may include editing the source code to extend, fix or replace parts of
PISM.

We provide both user-level (this manual) and developer-level documentation. Please see
source code browsers at |pism-docs| for the latter.

.. FIXME: mention debugging settings

- To use your (modified) version of PISM, you will need to follow the compilation from
  sources instructions in the :ref:`Installation Manual <sec-installation>`.
- We find it very useful to be able to check if a recent source code change broke
  something. PISM comes with "regression tests", which check if certain parts of PISM
  perform the way they should.\ [#]_

  Run "``make test``" in the build directory to run PISM's regression tests.

  Note, though, that while a test failure usually means that the new code needs more work,
  passing all the tests does not guarantee that everything works as it should. We are
  constantly adding new tests, but so far only a subset of PISM's functionality can be
  tested automatically.
- We strongly recommend using a version control system to manage code changes. Not only is
  it safer than the alternative, it is also more efficient.

.. rubric:: Footnotes

.. [#] This automates running verification tests described in section :ref:`sec-verif`,
       for example.
