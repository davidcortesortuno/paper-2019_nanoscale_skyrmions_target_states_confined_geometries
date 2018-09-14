#!/bin/bash
# Use a space in --B to avoid argparse to think that a negative
# number is another argument
for (( B = -2000; B < 2001; B += 100 )); do
    for (( R = 5; R < 6; R += 1 )); do
        SIMNAME="2Dhex_hexagon_R${R}nm_PdFe-Ir_1-tgt-st-up_B"
        echo "B=${B}e-3 T" "R=${R} nm"
        python ../hexagonal_fidimag.py \
                "${SIMNAME}${B}mT" \
                --hexagon ${R} 0.2715 \
                --D 1.557 \
                --B " ${B}e-3" \
                --J 5.881 \
                --mu_s 3 \
                --k_u 0.406 \
                --initial_state_skyrmion_up $(echo "$R - 1.0" | bc) 2 1 \
                --no_precession \
                --stopping_dmdt 0.000001 \
                --save_initial_vtk
    done
done
