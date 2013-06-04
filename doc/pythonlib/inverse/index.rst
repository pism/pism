==========================
PISM Inverse Problems
==========================

PISM contains algorithms for solving two closely related inverse problems:

  1. Given observed SSA velocity observations :math:`U_{\mathrm{obs}}`, determine effective yield stresses :math:`\tau_c`.
  2. Given observed SSA velocity observations :math:`U_{\mathrm{obs}}`, determine vertically averaged hardnesses :math:`B`.

The following sections contain an overview of the general inversion
approaches taken in PISM, a detailed description of the forward and 
inverse problems, and documentation of the :file:`pismi.py` 
command-line tool for running these algorithms.

.. toctree::
  refresher
  ssa_forward
  ssa_inverse
  pismi
  design