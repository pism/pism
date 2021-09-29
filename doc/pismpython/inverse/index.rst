=====================
PISM Inverse Problems
=====================

PISM contains algorithms for solving two closely related inverse problems:

  1. Given observed SSA velocity observations :math:`\vU_{\mathrm{obs}}`, determine effective yield stresses :math:`\tau_c`.
  2. Given observed SSA velocity observations :math:`\vU_{\mathrm{obs}}`, determine vertically averaged hardnesses :math:`B`.

The following sections contain an overview of the general inversion
approaches taken in PISM, a detailed description of the forward and 
inverse problems, and documentation of the :file:`pismi.py` 
command-line tool for running these algorithms.

We use the following typesetting conventions for variable names:

  * Config variables 
    (i.e. those appearing in :file:`pism_config.nc` or its overrides) 
    appear as: :cfg:`ssa_dirichlet_bc`.
  * Command line flags setting PISM/PETSc options appear the same way,
    but start with a hyphen: :cfg:`-ssa_dirichlet_bc`
  * Variables in NC files corresponding to data defined on a grid appear
    as: :ncvar:`vel_bc`.

.. toctree::
  refresher
  pism_ssa
  ssa_forward
  ssa_inverse
  pismi
  listeners
  design
