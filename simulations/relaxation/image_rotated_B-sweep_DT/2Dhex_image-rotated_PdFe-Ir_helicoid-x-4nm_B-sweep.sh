#!/bin/bash
# The argument --B needs an extra space at the beginning to avoid
# problems with dashes in the argument
IMG_PATH="../../../Romming_data/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw.png"
MU_S="../../../mu_s/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw_mus3-muB.npy"
#
for (( B = 0; B < 3001; B += 100 )); do
    SIMNAME="2Dhex_image-rotated_PdFe-Ir_helicoid-x-4nm_B${B}mT"
    python ../hexagonal_fidimag.py \
        --mesh_from_image "$IMG_PATH" 0 21.03 0 17.79 \
        "$SIMNAME" \
        --D 1.557 \
        --B " ${B}e-3" \
        --J 5.881 \
        --mu_s $MU_S \
        --k_u 0.406 \
        --initial_state_helicoid "4" "x" \
        --no_precession \
        --stopping_dmdt 0.000001 \
        --save_initial_vtk
done
