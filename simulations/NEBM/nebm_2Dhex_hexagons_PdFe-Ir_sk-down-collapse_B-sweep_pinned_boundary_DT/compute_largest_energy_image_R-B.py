import numpy as np
import argparse

# -----------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Find largest energy point')

parser.add_argument('energy_ndt',
                    help='Path to the ndt file with the energies of the '
                    'images in the band')

parser.add_argument('--first_saddle_point',
                    help='specify this option to return the first saddle '
                    'point', action='store_true')

# Parser arguments
args = parser.parse_args()

# -----------------------------------------------------------------------------

energies = np.loadtxt(args.energy_ndt)
# Get energy from last step and discard step number
energies = energies[-1][1:]
# Set the energy scale with respect to the energy of the first state
energies = energies - energies[0]

if not args.first_saddle_point:
    # Assuming there is no other point with the same energy
    # Get the array and then the element
    image_index = np.where(energies == energies.max())[0][0]
else:
    """
    We look for a saddle point with the profile:

                     i+2
                      O
                    /  \
                  -o    \
                 / i+1   o  i+3
                o         \
            .../ i         o i+4
                            \_ ...

    where the two closest neighbouring images in every direction
    have smaller energy
    """
    for i, E in enumerate(energies[2:-2]):
        if energies[i + 2] > energies[i + 1] and \
           energies[i + 1] > energies[i] and \
           energies[i + 2] > energies[i + 3] and \
           energies[i + 3] > energies[i + 4]:

            image_index = i + 2
            break

print(image_index)
