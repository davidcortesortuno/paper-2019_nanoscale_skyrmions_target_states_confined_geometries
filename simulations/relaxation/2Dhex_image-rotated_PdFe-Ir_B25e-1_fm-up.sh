#!/bin/bash
IMG_PATH="../../Romming_data/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw.png"
SIMNAME="2Dhex_image-rotated_PdFe-Ir_B25e-1_fm-up"
python hexagonal_fidimag.py \
        --mesh_from_image "$IMG_PATH" 0 21.03 0 17.79 \
        "$SIMNAME" \
        --D 1.557 \
        --B "2.5" \
        --J 5.881 \
        --mu_s 3 \
        --k_u 0.406 \
        --initial_state_ferromagnetic_up \
        --no_precession \
        --stopping_dmdt 0.0001
