.. include:: ../global.txt

.. _sec-citing-pism:

Citing PISM
===========

We ask that you cite the appropriate references if you publish results of a study that
involved running PISM or analyzing model outputs from a previously published PISM
application.

Receiving citations for PISM is important to demonstrate the relevance of our work to our
funding agencies and is a matter of fairness to everyone who contributed to its
development.\ [#f1]_

Unless PISM developers are involved in the preparation of a publication at the usual
co-author level, we do not expect co-authorship on papers using PISM.

Publications mentioning PISM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Please cite the following references in publications that mention but do not use PISM.\
[#f2]_

.. code-block:: bibtex
   :name: pism-citations
   :caption: Recommended PISM citations

   @article {BuelerBrown2009,
     AUTHOR = {E. Bueler and J. Brown},
      TITLE = {Shallow shelf approximation as a "sliding law" in a
               thermodynamically coupled ice sheet model},
    JOURNAL = {J. Geophys. Res.},
     VOLUME = {114},
        DOI = {10.1029/2008JF001179},
       YEAR = {2009},
   }

   @article {Winkelmannetal2011,
     AUTHOR = {Winkelmann, R. and Martin, M. A. and Haseloff, M. and Albrecht, T.
               and Bueler, E. and Khroulev, C. and Levermann, A.},
      TITLE = {The {P}otsdam {P}arallel {I}ce {S}heet {M}odel ({PISM-PIK})
               {P}art 1: {M}odel description},
    JOURNAL = {The Cryosphere},
     VOLUME = {5},
       YEAR = {2011},
      PAGES = {715--726},
        DOI = {10.5194/tc-5-715-2011},
   }

Publications using PISM
%%%%%%%%%%%%%%%%%%%%%%%

If you use PISM or model outputs from an earlier PISM application in your publication then
we ask for :ref:`an acknowledgment of funding <sec-funding-acknowledgment>` and a
citation. Please cite :ref:`papers listed above <pism-citations>` as well as the code
itself. See :ref:`sec-code-availability` for details.

.. note::

   We maintain a list of PISM applications at `pism.io/publications
   <https://www.pism.io/publications>`_. Please send an e-mail to |pism-email| to let us
   know about your work!

See the `publishing guidelines <geodynamics-publishing-guidelines_>`_ from Computational
Infrastructure for Geodynamics (CIG) for more recommendations.

.. _sec-funding-acknowledgment:

Funding acknowledgment
++++++++++++++++++++++

Please include this statement to acknowledge PISM's funding:

    .. include:: ../funding.txt

.. _sec-code-availability:

Code availability
+++++++++++++++++

Please state the exact PISM version used in your study and cite the corresponding Zenodo
record; see :numref:`tab-pism-doi` or https://doi.org/10.5281/zenodo.1199019 to find the
correct DOI. You can then export citation information from a record's web page.

.. admonition:: Example statement

   We use PISM version 2.1, which is available under the GPL license from Zenodo
   [cite the appropriate Zenodo record].

Please use a DOI instead of an URL pointing to |pism-website| or |pism-github-url|
whenever possible.

Please cite |pism-website| when talking about code that will be available in *future*
PISM versions.

You can cite all PISM versions by using the DOI `10.5281/zenodo.1199019
<https://doi.org/10.5281/zenodo.1199019>`_. This DOI represents all versions and will
always resolve to the latest one.

.. csv-table:: DOIs for recent PISM versions
   :header: Version, DOI
   :name: tab-pism-doi

   2.1,   `10.5281/zenodo.10202029 <https://doi.org/10.5281/zenodo.10202029>`_
   2.0.7, `10.5281/zenodo.10183040 <https://doi.org/10.5281/zenodo.10183040>`_
   2.0.6, `10.5281/zenodo.7563365  <https://doi.org/10.5281/zenodo.7563365>`_
   2.0.5, `10.5281/zenodo.7199611  <https://doi.org/10.5281/zenodo.7199611>`_
   2.0.4, `10.5281/zenodo.6584915  <https://doi.org/10.5281/zenodo.6584915>`_
   2.0.3, `10.5281/zenodo.6537481  <https://doi.org/10.5281/zenodo.6537481>`_
   2.0.2, `10.5281/zenodo.6001196  <https://doi.org/10.5281/zenodo.6001196>`_
   2.0.1, `10.5281/zenodo.5966037  <https://doi.org/10.5281/zenodo.5966037>`_
   2.0,   `10.5281/zenodo.5758039  <https://doi.org/10.5281/zenodo.5758039>`_

Studies using customized versions of PISM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We strongly encourage you to archive your code using Zenodo_ or a similar service and
obtain a DOI for the resulting record. Please make sure to cite this record **and** the
PISM version your work is based on.

.. admonition:: Example statement

   We use a modified version of PISM version 2.1 [cite PISM 2.1] which is available under
   the GPL license from ``<Research data repository name>`` [citation for the
   corresponding record].

Uploading the code to a public repository (e.g. on GitHub_) and providing a URL in a
publication is not sufficient: repository contents may change, a repository may be
deleted, or a company hosting it may stop providing this service.


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

.. [#f2] This is appropriate when providing a citation for a method first employed in
         PISM, e.g. "The ice velocity is calculated as the sum of the SIA and SSA
         contributions [PISM citation]."
