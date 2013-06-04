.. _SSAInverse:

SSA Inverse Problem
===================

PISM supports two inverse problems for the SSA.

  1. Given observed SSA velocity observations :math:`U_{\mathrm{obs}}`, determine effective yield stress :math:`\tau_c`.
  2. Given observed SSA velocity observations :math:`U_{\mathrm{obs}}`, determine vertically averaged hardness :math:`B`.
  
These are solved using ``pismi -inv_ssa tauc`` and 
``pismi -inv_ssa hardav`` respectively.
