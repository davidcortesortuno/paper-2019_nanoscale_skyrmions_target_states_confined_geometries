# import numpy as np
import textwrap
import subprocess
from subprocess import Popen
import sys


initial_states = [[r'initial_state_skyrmion_down $(echo "$R * 0.1 / 4.0" | bc)', 'SKD'],
                  [r'initial_state_skyrmion_up $(echo "$R * 0.1 - 1.0" | bc) 2 1', 'TGTU'],
                  [r'initial_state_skyrmion_down $(echo "$R * 0.1 - 1.0" | bc) 3 1', '3PID'],
                  [r'initial_state_skyrmion_up $(echo "$R * 0.1 - 1.0" | bc) 4 1', '4PIU'],
                  [r'initial_state_ferromagnetic_up', 'FMU'],
                  [r'initial_state_helicoid "2" "x"', 'HX2'],
                  [r'initial_state_helicoid "4" "x"', 'HX4'],
                  [r'initial_state_helicoid "6" "x"', 'HX6'],
                  [r'initial_state_helicoid "2" "y"', 'HY2'],
                  [r'initial_state_helicoid "4" "y"', 'HY4'],
                  [r'initial_state_helicoid "6" "y"', 'HY6'],
                  [r'initial_state_random 42', 'R42'],
                  # [r'initial_state_random 4242', 'R4242'],
                  [r'initial_state_random 24', 'R24'],
                  # [r'initial_state_random 2424', 'R2424'],
                  [r'initial_state_random 99', 'R99'],
                  # [r'initial_state_random 9999', 'R9999'],
                  [r'initial_state_multiple_dots_down_REL 0.05 0.25 0.0 -0.25 0.0', 'MD2'],
                  [r'initial_state_multiple_dots_down_REL 0.05 0.25 -0.15 -0.25 -0.15 0.0 0.2', 'MD3'],
                  [r'initial_state_multiple_dots_down_REL 0.05 0.25 -0.25 0.25 0.25 -0.25 0.25 -0.25 -0.25', 'MD4'],
                  [r'initial_state_multiple_dots_down_REL 0.05 '
                   r'0.3 -0.3 0.3 0.3 -0.3 0.3 -0.3 -0.3 0.0 0.0', 'MD5'],
                  [r'initial_state_multiple_dots_down_REL 0.05 '
                   r'0.3 -0.3 0.3 0.3 -0.3 0.3 -0.3 -0.3 0.0 0.6 0.0 0.0', 'MD6'],
                  [r'initial_state_multiple_dots_down_REL 0.05 '
                   r'0.3 -0.3 0.3 0.3 -0.3 0.3 -0.3 -0.3 0.0 0.6 0.0 -0.6 0.0 0.0', 'MD7'],
                  [r'initial_state_n_dots_down_helicoid_REL 0.4 0.5 0.05 0.05', 'MD1H'],
                  [r'initial_state_n_dots_down_helicoid_REL 0.4 0.7 0.4 0.3 0.05 0.05', 'MD2H'],
                  [r'initial_state_skyrmion_down $(echo "$R * 0.1 - 1.0" | bc) 5 1', '5PID'],
                  ]
# 23 Initial States

IS = initial_states[int(sys.argv[1])]
job = ("""\
#!/bin/bash
for (( B = 0; B < 2001; B += 50 )); do
    for (( R = 40; R < 201; R += 5 )); do
        SIMNAME="PdFe-Ir_hexagon_R${{R}}e-1nm_B${{B}}mT_{0}"
        echo "B=${{B}}e-3 T" "R=${{R}}e-1 nm"
        python ../hexagonal_fidimag.py \\
            "${{SIMNAME}}" \\
            --hexagon ${{R}}e-1 0.2715 \\
            --D 1.557 \\
            --B " ${{B}}e-3" \\
            --J 5.881 \\
            --mu_s 3 \\
            --k_u 0.406 \\
            --{1} \\
            --no_precession \\
            --stopping_dmdt 0.000001 \\
            --pin_boundaries 0 0 1 \\
            --save_energy \\
            --save_Q
done
done
""".format(IS[1], IS[0])
)
# print(job)
f = open('{}.sh'.format(IS[1]), 'w')
f.write(job)
f.close()
subprocess.call('bash {}.sh'.format(IS[1]), shell=True)
