.. include:: global.txt

.. _sec-publishing:

Publishing
==========

If you use PISM or model outputs from an earlier PISM application in your publication then
we ask for an acknowledgment of funding (see :ref:`sec-funding-acknowledgment`) and a
citation. Please see :ref:`sec-citing-pism` for details.

Receiving citations for PISM is important to demonstrate the relevance of our work to our
funding agencies and is a matter of fairness to everyone who contributed to its
development.\ [#f1]_

We do not expect co-authorship on papers using PISM unless PISM developers are involved in
the preparation of a publication at the usual co-author level.

.. note::

   We maintain a list of PISM applications at `pism.io/publications
   <https://www.pism.io/publications>`_. Please send an e-mail to |pism-email| to let us
   know about your work!

See the `publishing guidelines <geodynamics-publishing-guidelines_>`_ from Computational
Infrastructure for Geodynamics (CIG) for more recommendations.

.. _sec-funding-acknowledgment:

Funding acknowledgment
%%%%%%%%%%%%%%%%%%%%%%

Please include this statement to acknowledge PISM's funding:

    .. include:: funding.txt

Open research
%%%%%%%%%%%%%

We encourage you to make code and data used in your study publicly available to simplify
reproducing your results and provide sufficient information to reviewers to evaluate the
publication.

.. _sec-code-availability:

Code availability
+++++++++++++++++

If you used a released PISM version *without modifications* it is sufficient to :ref:`cite
the specific PISM version <sec-citing-pism>` you used.

If you used a *development* version of PISM (even if you did not make any modifications),
please archive and cite it as if it is a :ref:`customized PISM version
<sec-citing-customized-pism>`.

.. admonition:: Example statement

   We use PISM version 2.1, which is available under the GPL license from Zenodo
   [cite the appropriate Zenodo record].

.. _sec-citing-customized-pism:

Studies using customized versions of PISM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We strongly encourage you to archive your code using a research data repository (e.g.
Zenodo_) and obtain a DOI for the resulting record. Please make sure to cite this record
**and** the specific PISM version your work is based on.

.. admonition:: Example statement

   We use a modified version of PISM version 2.1 [cite PISM 2.1] which is available under
   the GPL license from ``<Research data repository name>`` [citation for the
   corresponding record].

Uploading the code to a public repository (e.g. on GitHub_) and providing a URL in a
publication is not sufficient: repository contents may change, a repository may be
deleted, or a company hosting it may stop providing this service.

.. note::

   Please consider contributing your modifications to PISM's repository! See
   :ref:`sec-contributing` for details.


.. _sec-data-availability:

Data availability
+++++++++++++++++

We recommend archiving

- all input files necessary to reproduce your research,
- all pre-processing scripts that were used to prepare model inputs,
- shell scripts and configuration files needed to run PISM,
- PISM outputs resulting from simulations you performed,
- any output files that can be used to verify that the reproduction of results was
  successful,
- post-processing, analysis and plotting scripts used in preparation of a publication,
- a list of required pre-processing, post-processing and analysis tools, including their
  versions

and then citing the corresponding record using its DOI.

Sometimes an input data set is not publicly available and cannot be archived this way. If
possible, please explain how to request access to it (e.g. "please contact *Author X et
al* to request access").

If model outputs are prohibitively large (e.g. if you ran an ensemble of thousands of
high-resolution simulations), please select and archive a subset of outputs that can be
used as a "standard." A person attempting to reproduce your study can then attempt to
reproduce these outputs to verify their model setup.

.. admonition:: Example statement

   The PISM code, model inputs and outputs, the scripts used to run simulations, analyze
   results and create figures are available at ``https://doi.org/<insert DOI>`` (cite the
   DOI you obtained).

.. caution::

   Each PISM output file includes the PISM version, versions of major PISM dependencies
   and values of all configuration parameters that were used. Please ensure that this
   information is not stripped from archived files during post-processing.

.. rubric:: Footnotes

.. [#f1] See :ref:`sec-pism-authors` for the list of authors and their contributions.
