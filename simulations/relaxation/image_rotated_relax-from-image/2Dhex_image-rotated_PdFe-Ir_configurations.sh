#!/bin/bash
# Script to simulate a hexagonal 2D system of FePd over Ir
# atoms. Parameters are taken from Romming et al
# publications
#
BS=("0" "1" "1.25" "1.61" "2" "2.5" "2.5" "3" "4" "1.60" "1.60" "-0.5" "-0.5" "-1.25" "-1.60")
#
for i in `seq 0 14`
do
    NUM=$(printf "%03d" ${i})
    IMG_PATH="../../../Romming_data/Pd_Fe_Ir111_configurations/gimp/png/${NUM}.png"
    SIMNAME="2Dhex_image-rotated_PdFe-Ir_configuration-${NUM}"
    MU_S="../../../mu_s/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw_mus3-muB.npy"
    python ../hexagonal_fidimag.py \
            "$SIMNAME" \
            --mesh_from_image "$IMG_PATH" 0 21.03 0 17.79 \
            --D 1.557 \
            --B " ${BS[${i}]}" \
            --J 5.881 \
            --mu_s 3 \
            --k_u 0.406 \
            --initial_state_from_image \
            --no_precession \
            --stopping_dmdt 0.00001
done
