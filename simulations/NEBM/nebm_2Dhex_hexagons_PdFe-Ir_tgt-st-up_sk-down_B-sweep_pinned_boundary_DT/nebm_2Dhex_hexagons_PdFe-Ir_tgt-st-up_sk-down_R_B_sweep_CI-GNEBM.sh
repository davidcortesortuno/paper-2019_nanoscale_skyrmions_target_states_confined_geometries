#!/bin/bash 
#
mkdir -p ndts

simulation () {
    R=$1
    B=$2
    SIMNAME="nebm_2Dhex_hexagon_R${R}nm_PdFe-Ir_tgt-st-up_sk_down_B${B}mT_CI-GNEBM"
    ORIG_SIMNAME="nebm_2Dhex_hexagon_R${R}nm_PdFe-Ir_tgt-st-up_sk_down_B${B}mT"
    CLIMB_IMAGE=$(python compute_largest_energy_image_R-B.py ndts/${ORIG_SIMNAME}_energy.ndt)
    python ../hexagonal_neb_fidimag.py \
        --sim_name $SIMNAME \
        --hexagon ${R} 0.2715 \
        --D 1.557 \
        --B " ${B}e-3" \
        --J 5.881 \
        --mu_s 3 \
        --k_u 0.406 \
        --images_files npys/${ORIG_SIMNAME}_LAST \
        --nebm_k 1e4 \
        --nebm_steps 2000 --save_vtk 2000 --save_npy 2000 \
        --stopping_dYdt 0.00001 \
        --save_interpolation 500 \
        --climbing_image ${CLIMB_IMAGE} \
        --pin_boundaries 0 0 1

    mv ${SIMNAME}_energy.ndt ndts/
    mv ${SIMNAME}_dYs.ndt ndts/
}

N=4
for R in 6 8 10 12 14; do
    for (( B = 0; B < 1201; B += 100 )); do
        if ((R < 10  || B > 800)); then
            ((i=i%N))
            ((i++==0)) && wait
            simulation $R $B &!
        fi
    done
done
