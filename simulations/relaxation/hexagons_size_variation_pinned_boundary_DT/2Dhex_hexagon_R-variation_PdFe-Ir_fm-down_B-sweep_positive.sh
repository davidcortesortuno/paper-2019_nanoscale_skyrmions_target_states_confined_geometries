#!/bin/bash
# Use a space in --B to avoid argparse to think that a negative
# number is another argument
for (( B = -2000; B < 1; B += 100 )); do
    for (( R = 22; R < 41; R += 2 )); do
        SIMNAME="2Dhex_hexagon_R${R}nm_PdFe-Ir_fm-down_B"
        echo "B=${B}e-3 T" "R=${R} nm"
        python ../hexagonal_fidimag.py \
                "${SIMNAME}${B}mT" \
                --hexagon ${R} 0.2715 \
                --D 1.557 \
                --B " ${B}e-3" \
                --J 5.881 \
                --mu_s 3 \
                --k_u 0.406 \
                --initial_state_ferromagnetic_down \
                --no_precession \
                --stopping_dmdt 0.000001 \
                --pin_boundaries 0 0 -1
    done
done
