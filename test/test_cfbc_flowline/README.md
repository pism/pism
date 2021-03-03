### Verification test for modified driving stress at calving front using van-der-Veen solution

- Copyright (C) 2021 Ronja Reese, Torsten Albrecht and PISM authors
- contact: reese@pik-potsdam.de and albrecht@pik-potsdam.de

This is analogous to `examples/marine/flowline` or `test/test_shelf/exactV.py`
Simply run 'bash run_iceshelf_test.sh'


The plots show, that altering the ice thickness in the last grid cell at the calving front only (e.g. local melting) affects the CFBC and hence the entire SSA solution. The driving stress scheme has been corrected accordingly. Without perturbation both schemes agree well with the unbuttressed steady state solution like van der Veen.

More details will follow in **Reese et al., in prep.**
