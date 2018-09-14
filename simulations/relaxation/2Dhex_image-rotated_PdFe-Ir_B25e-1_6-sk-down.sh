
#!/bin/bash
IMG_PATH="../../Romming_data/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw.png"
MU_S="../../mu_s/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw_mus3-muB.npy"
SIMNAME="2Dhex_image-rotated_PdFe-Ir_B25e-1_6-sk-down"
python hexagonal_fidimag.py \
        --mesh_from_image "$IMG_PATH" 0 21.03 0 17.79 \
        "$SIMNAME" \
        --D 1.557 \
        --B "-2.5" \
        --J 5.881 \
        --mu_s $MU_S \
        --k_u 0.406 \
        --initial_state_multiple_sks_down 2 7.5 5.0 \
                                          7.5 11.5 \
                                          14.5 13.0 \
                                          14.5 5.0 \
                                          11 9.5 \
                                          11 2.5 \
        --no_precession \
        --stopping_dmdt 0.000001 \
        --save_initial_vtk
