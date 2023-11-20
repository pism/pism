### Scripts using circular shelf setups

These scripts use simplified circular shelf setups to test

* the parameterization of the sub-grid scale motion of the calving front (`run_noshelf.sh`, `run_shelfonly.sh`, `run_withshelf.sh`)
* the "calving at a threshold thickness" mechanism (`test_calving_at_thickness.sh`)
* the "eigencalving" mechanism (`test_eigencalving.sh`)
* the "iceberg removal" machanism (`test_iceberg_removal.sh`; this one is not circular, though)

The script `circular_dirichlet.py` generates input files with the prescribed ice flux at the grounding line and is used both here and elsewhere (e.g. in `examples/marine/melange`).
