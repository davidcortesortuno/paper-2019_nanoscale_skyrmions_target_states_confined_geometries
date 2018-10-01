#!/bin/bash
# The argument --B needs an extra space at the beginning to avoid
# problems with dashes in the argument
IMG_PATH="../../../Romming_data/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw.png"
MU_S="../../../mu_s/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw_mus3-muB.npy"

mkdir -p ndts
for (( B = 0; B < 501; B += 100 )); do
    echo "Field: ${B} mT ** Climbing Image: ${CI}"
    #
    ORIG_SIMNAME="nebm_2Dhex_image-rotated_PdFe-Ir_3PI-vortex_2PI-vortex_B${B}mT"
    SIMNAME="${ORIG_SIMNAME}_CI-GNEBM"

    CLIMB_IMAGE=$(python compute_largest_energy_image.py ndts/${ORIG_SIMNAME}_energy.ndt)
    #
    python ../hexagonal_neb_fidimag.py \
        --sim_name $SIMNAME \
        --mesh_from_image "$IMG_PATH" 0 21.03 0 17.79 \
        --D 1.557 \
        --B " ${B}e-3" \
        --J 5.881 \
        --mu_s $MU_S \
        --k_u 0.406 \
        --images_files npys/${ORIG_SIMNAME}_LAST \
        --nebm_k 1e4 \
        --nebm_steps 2000 --save_vtk 2000 --save_npy 2000 \
        --stopping_dYdt 0.00001 \
        --climbing_image ${CLIMB_IMAGE} \
        --save_interpolation 500
    mv ${SIMNAME}_energy.ndt ndts/
    mv ${SIMNAME}_dYs.ndt ndts/
done
