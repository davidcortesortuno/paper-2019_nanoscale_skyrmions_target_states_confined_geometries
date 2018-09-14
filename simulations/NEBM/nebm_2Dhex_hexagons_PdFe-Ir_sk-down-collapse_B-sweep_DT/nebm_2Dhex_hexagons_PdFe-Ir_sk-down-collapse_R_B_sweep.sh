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
        --interpolation "../../relaxation/hexagons_size_variation_DT/npys/2Dhex_hexagon_R${R}nm_PdFe-Ir_1-sk-down_B${B}mT_npys" \
                         26 \
                         "../../relaxation/hexagons_size_variation_DT/npys/2Dhex_hexagon_R${R}nm_PdFe-Ir_fm-up_B${B}mT_npys" \
        --nebm_k 1e4 \
        --nebm_steps 2000 --save_vtk 2000 --save_npy 2000 \
        --stopping_dYdt 0.0001 \
        --save_interpolation 500
    mv ${SIMNAME}_energy.ndt ndts/
    mv ${SIMNAME}_dYs.ndt ndts/
}

N=4
for R in 3 4 6 8 10; do
    for (( B = 0; B < 2001; B += 100 )); do
        ((i=i%N))
        ((i++==0)) && wait
        simulation $R $B &!
    done
done
