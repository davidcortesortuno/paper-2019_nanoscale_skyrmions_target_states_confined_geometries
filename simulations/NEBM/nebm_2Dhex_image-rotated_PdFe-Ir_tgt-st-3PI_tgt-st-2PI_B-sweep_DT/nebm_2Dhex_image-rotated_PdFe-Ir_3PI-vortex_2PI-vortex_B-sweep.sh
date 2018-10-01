#!/bin/bash
# The argument --B needs an extra space at the beginning to avoid
# problems with dashes in the argument
IMG_PATH="../../../Romming_data/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw.png"
MU_S="../../../mu_s/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw_mus3-muB.npy"

mkdir -p ndts
for (( B = 0; B < 501; B += 100 )); do
    SIMNAME="nebm_2Dhex_image-rotated_PdFe-Ir_3PI-vortex_2PI-vortex_B${B}mT"
    #
    python ../hexagonal_neb_fidimag.py \
        --sim_name $SIMNAME \
        --mesh_from_image "$IMG_PATH" 0 21.03 0 17.79 \
        --D 1.557 \
        --B " ${B}e-3" \
        --J 5.881 \
        --mu_s ${MU_S} \
        --k_u 0.406 \
        --interpolation "../../relaxation/image_rotated_B-sweep_DT/npys/2Dhex_image-rotated_PdFe-Ir_tgt-st-up-3PI_B${B}mT_npys/" \
                        28 \
                        "../../relaxation/image_rotated_B-sweep_DT/npys/2Dhex_image-rotated_PdFe-Ir_tgt-st-down_B${B}mT_npys/" \
        --nebm_k 1e4 \
        --nebm_steps 2000 --save_vtk 2000 --save_npy 2000 \
        --stopping_dYdt 0.0001 \
        --save_interpolation 500
    #
    mv ${SIMNAME}_energy.ndt ndts/
    mv ${SIMNAME}_dYs.ndt ndts/
done
