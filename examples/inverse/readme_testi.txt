compare

python invert_testi_tik.py -tao_fatol 1e-6 -tao_frtol 1e-10 -tao_monitor -pseudo_plastic -pseudo_plastic_q 0. -inv_ssa_tauc_param trunc -inv_ssa_cL2 0 -inv_ssa_cH1 1 -eta 1e-5 -tao_gttol 1e-4 -tauc_guess_scale .3
 
python inv_testi_basic.py -eta 1e-5 -tao_fatol 1e-6 -tao_frtol 1e-10 -tao_monitor -is_dimensionless 0 -inv_ssa_tauc_param exp -inv_ssa_cL2 0 -inv_ssa_cH1 1 -tau_gttol 1e-4 -tauc_guess_scale .3

which should both result in 20% RMS error
