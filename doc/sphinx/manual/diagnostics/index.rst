.. include:: ../../global.txt

.. To update this list, run PISM with the option -list_diagnostics_json, making sure that
   you enabled all the sub-models that provide diagnostics (scalar and
   spatially-variable).

   At the time of writing this includes

   - Evolving geometry (without -no_mass)
   - The bedrock thermal layer model (-Mbz X for X > 1)
   - An energy-conservation model (-energy enthalpy or -energy cold)
   - A bed deformation model (-bed_def iso or -bed_def lc)
   - The PDD surface model (-surface pdd)
   - An atmosphere model using a cosine yearly cycle (e.g. -atmosphere searise_greenland)
   - The hybric stress balance model (-stress_balance ssa+sia -ssa_method fd)
   - The Mohr-Coulomb basal yield stress model (-yield_stress mohr_coulomb)
   - The "routing" subglacial hydrology model (-hydrology routing)
   - All supported calving mechanisms (-calving
     eigen_calving,thickness_calving,frontal_melt,vonmises_calving)

   The "distributed" hydrology model does add more diagnostics, but it is not supported at
   this point (and not documented in the User's Manual).

.. _sec-diagnostics-list:

Diagnostic quantities
=====================

The availability of a diagnostic depends on modeling choices. For example, the bed uplift
rate :var:`dbdt` is available only if a bed deformation model is selected.

Some scalar diagnostics come in two versions: the one with the suffix ``_glacierized`` and
the one without. Here the former account for the ice
:config:`output.ice_free_thickness_standard` meters or thicker (10 meters by default) and
the latter include all ice regardless of the thickness. "Glacierized" versions were added
to make it easier to analyze changes in glacier volumes and areas and exclude changes in
the seasonal snow cover.

.. contents::

.. include:: diagnostics-list.txt
