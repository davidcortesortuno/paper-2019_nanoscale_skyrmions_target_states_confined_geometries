#!/bin/bash 
#
mkdir -p ndts

simulation () {
    R=$1
    B=$2
    SIMNAME="nebm_2Dhex_hexagon_R${R}nm_PdFe-Ir_sk_down-fm_up_B${B}mT"
    python ../hexagonal_neb_fidimag.py \
        --sim_name $SIMNAME \
        --hexagon ${R} 0.2715 \
        --D 1.557 \
        --B " ${B}e-3" \
        --J 5.881 \
        --mu_s 3 \
        --k_u 0.406 \
        --images_files initial_states/npys/R${R}nm_B${B}mT \
        --nebm_k 1e4 \
        --nebm_steps 2000 --save_vtk 2000 --save_npy 2000 \
        --stopping_dYdt 0.000001 \
        --save_interpolation 500
    mv ${SIMNAME}_energy.ndt ndts/
    mv ${SIMNAME}_dYs.ndt ndts/
}

N=4
for R in 10 8 6 4 3; do
    for (( B = 700; B < 800; B += 100 )); do
        if ((R != 10 || B >= 700)); then
            ((i=i%N))
            ((i++==0)) && wait
            simulation $R $B &!
        fi
    done
done
