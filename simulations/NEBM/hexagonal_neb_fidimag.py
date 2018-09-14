from __future__ import print_function
import argparse

# ARGUMENTS
parser = argparse.ArgumentParser(description='NEB method for 2D films with '
                                 'interfacial DMI')
method = parser.add_mutually_exclusive_group(required=True)

mesh = parser.add_mutually_exclusive_group(required=True)

mesh.add_argument('--hexagonal_mesh',
                  help='Generates a mesh using a hexagonal '
                  'lattice. Arguments are: [number_spins_x number_spins_y '
                  'lattice_constant]. The spins alignment can be changed '
                  'with the corresponding optional argument',
                  type=float, nargs=3)

mesh.add_argument('--truncated_triangle',
                  help='Generates a truncated triangle geometry using a '
                  ' hexagonal lattice. Arguments are: '
                  ' [triangle_side_length offsets lattice_const] '
                  'The offsets can be a float or a 3 element list',
                  type=float, nargs='+')

mesh.add_argument('--hexagon',
                  help='Generates a hexagon geometry using a '
                  ' hexagonal lattice. Arguments are: '
                  ' [hexagon_circumradius lattice_const] ',
                  type=float, nargs=2)

mesh.add_argument('--mesh_from_image',
                  help='Generates the mesh using an '
                  'image as a reference for the position of the magnetic '
                  'moments in a hexagonal lattice. The arguments are: '
                  '[image_path x_min x_max y_min y_max]. The last four are '
                  'the limits of the image given in nm. This option uses the'
                  'sim_from_image library',
                  nargs=5
                  )

# NEB Simulation --------------------------------------------------------------

parser.add_argument('--nebm_steps', help='Number of steps for the '
                    'Cartesian or Spherical NEB method')

parser.add_argument('--save_vtk', help='Number specifying that the vtk files'
                    ' are going to be saved every *save_vtk* number of steps',
                    type=int)

parser.add_argument('--save_npy', help='Number specifying that the npy files '
                    'are going to be saved every *save_vtk* number of steps',
                    type=int)

parser.add_argument('--sim_name',
                    help='Simulation name')

parser.add_argument('--unit_length', help='Mesh unit length',
                    type=float, default=1.0)

parser.add_argument('--alpha', help='Damping constant value',
                    type=float)

parser.add_argument('--PBC_2D',
                    help='Two dimensional boundary condition',
                    action='store_true')

parser.add_argument('--nebm_k', help='Spring constant for the NEB method'
                    '. Default: k=1e4',
                    default='1e4', type=float)

parser.add_argument('--stopping_dYdt', help='Stopping dYdt for the NEB method'
                    '. Default: 1e-2',
                    default='1e-2', type=float)

parser.add_argument('--tols', help='Tolerances for the integrator as '
                    'a pair: rtol atol',
                    type=float, nargs=2)

# Material  -------------------------------------------------------------------

parser.add_argument('--D', help='DMI constant in units of meV',
                    type=float, default=1)

parser.add_argument('--J', help='Exchange constant in units of meV',
                    type=float, default=1)

parser.add_argument('--mu_s', help='Magnetisation in units of mu_B',
                    default=2)

parser.add_argument('--k_u', help='Anisotropy constant in units of meV',
                    type=float)

parser.add_argument('--B', help='External magnetic field perpendicular to the'
                    ' square plane (z direction), in Tesla',
                    type=float)

parser.add_argument('--Demag', help='Add this option to use dipolar '
                    'interactions (TESTING)',
                    action='store_true')

parser.add_argument('--alignment', help='Hexagonal lattice alignment '
                    '(square or diagonal)',
                    default='diagonal')

parser.add_argument('--pin_boundaries', help='Set this option to fix the'
                    'spin directions at the boundaries. It must be specified '
                    'the three directions of the spins at this region',
                    nargs=3, type=float)

# NEB Method ------------------------------------------------------------------

parser.add_argument('--climbing_images',
                    help='An optional list of integers with values ranging '
                    'from 1 to the total '
                    ' number of images minus two (the extremes do not count)'
                    ' in order to use those images of the band with the'
                    ' Climbing Image NEB method. These images must be the'
                    ' maximum energy states and will not be affected by'
                    ' the Spring Force. This helps to have a better'
                    ' resolution around saddle points',
                    # metavar=('INDEX'),
                    type=int, nargs='+'
                    )

# Define a method: interpolation or images_files for the initial state
# which are mutually exclusive
method.add_argument('--interpolation',
                    help='Use a series of images for the initial state, '
                    'interpolating the magnetisation vectors components '
                    'from the first argument until the last argument.'
                    'The parameters are: '
                    'file1 number_interps file2 number_interps file3 ...'
                    'Where number_interps is the number of interpolations'
                    'between the i and (i+1)th states of the Energy band.'
                    ' Thus the number of arguments is always odd.',
                    nargs='+',
                    # metavar=('FILE1', 'N_INTERPS', 'FILE2', 'etc')
                    )

method.add_argument('--interpolation_rotation',
                    help='Use Rodrigues formula for interpolation ',
                    nargs='+'
                    )

method.add_argument('--images_files',
                    help='Use all the *npy* files from the specified folder '
                    'path for the initial states of the Energy Band. File '
                    'names must be sorted as *image_{}.pny* with {} as the '
                    'number in the sequence of images in the Energy Band, '
                    'and include the initial and final states',
                    # metavar=('NPY_FILES_PATH')
                    )

parser.add_argument('--save_interpolation', help='Save a dat file with a '
                    'cubic interpolation of the last NEBM step. Argument is '
                    'the number of data points',
                    type=int, default=100)

# -----------------------------------------------------------------------------

# Parser arguments
args = parser.parse_args()

import sys
import numpy as np

from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import UniformExchange, DMI, Anisotropy, DemagHexagonal
from fidimag.atomistic import Zeeman

# Import physical constants from fidimag
import fidimag.common.constant as const

from fidimag.common.nebm_geodesic import NEBM_Geodesic

import re
import os, glob
import shutil

this_path = os.path.abspath(__file__)
dir_this_path = os.path.dirname(this_path)

sys.path.append(dir_this_path + '/../../')
import sim_from_image as sfi

sys.path.append(dir_this_path + '/../../mesh_geometries')
from sim_geometries import TruncatedTriangleSim, HexagonSim

mu0 = 4 * np.pi * 1e-7


# -----------------------------------------------------------------------------
# Mesh ------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if args.hexagonal_mesh:
    if not args.PBC_2D:
        mesh = HexagonalMesh(args.hexagonal_mesh[2] * 0.5,
                             int(args.hexagonal_mesh[0]),
                             int(args.hexagonal_mesh[1]),
                             alignment=args.alignment,
                             unit_length=args.unit_length
                             )
        sim = Sim(mesh, name=args.sim_name)

    else:
        mesh = HexagonalMesh(args.hex_a, args.hex_nx, args.hex_ny,
                             periodicity=(True, True),
                             alignment=args.alignment,
                             unit_length=args.unit_length
                             )

    # Initiate Fidimag simulation ---------------------------------------------
    sim = Sim(mesh, name=args.sim_name)

elif args.mesh_from_image:
    sim_from_image = sfi.sim_from_image(args.mesh_from_image[0],
                                        image_range=[float(args.mesh_from_image[1]),
                                                     float(args.mesh_from_image[2]),
                                                     float(args.mesh_from_image[3]),
                                                     float(args.mesh_from_image[4])
                                                     ],
                                        sim_name=args.sim_name
                                        )

    if args.mu_s.endswith('.npy'):
        sim_from_image.generate_magnetic_moments(load_file=str(args.mu_s))
    else:
        sim_from_image.generate_magnetic_moments(mu_s=(float(args.mu_s) * const.mu_B))

    sim = sim_from_image.sim

elif args.truncated_triangle:
    if len(args.truncated_triangle) == 3:
        sim_triangle = TruncatedTriangleSim(
            args.truncated_triangle[0],  # L
            args.truncated_triangle[1],  # offset
            args.truncated_triangle[2],  # a
            float(args.mu_s),            # mu_s
            name=args.sim_name
            )
    elif len(args.truncated_triangle) == 5:
        sim_triangle = TruncatedTriangleSim(
            args.truncated_triangle[0],    # L
            [float(offs) for offs in args.truncated_triangle[1:4]],  # offsets
            args.truncated_triangle[4],    # a
            float(args.mu_s),              # mu_s
            name=args.sim_name
            )

    sim = sim_triangle.sim


elif args.hexagon:
    sim_hexagon = HexagonSim(args.hexagon[0],   # R
                             args.hexagon[1],   # a
                             float(args.mu_s),  # mu_s
                             name=args.sim_name
                             )
    sim = sim_hexagon.sim

# sim.set_tols(rtol=1e-10, atol=1e-14)
if args.alpha:
    sim.driver.alpha = args.alpha
# sim.gamma = 2.211e5

# Material parameters ---------------------------------------------------------

if args.hexagonal_mesh:
    sim.mu_s = args.mu_s * const.mu_B

exch = UniformExchange(args.J * const.meV)
sim.add(exch)

dmi = DMI(D=(args.D * const.meV), dmi_type='interfacial')
sim.add(dmi)

if args.B:
    zeeman = Zeeman((0, 0, args.B))
    sim.add(zeeman, save_field=True)

if args.k_u:
    # Uniaxial anisotropy along + z-axis
    sim.add(Anisotropy(args.k_u * const.meV, axis=[0, 0, 1]))

if args.Demag:
    sim.add(DemagHexagonal())

# -------------------------------------------------------------------------

# We will correct the spin directions according to the specified argument,
# in case the spins at the boundaries are pinned
if args.pin_boundaries:
    boundary_spins = np.logical_and(np.any(sim.mesh.neighbours < 0, axis=1),
                                    sim.mu_s != 0)

    new_m = np.copy(sim.spin.reshape(-1, 3))
    new_m[boundary_spins] = np.array([args.pin_boundaries[0],
                                      args.pin_boundaries[1],
                                      args.pin_boundaries[2]
                                      ])
    sim.set_m(new_m.reshape(-1,))

    # Now we pin the spins:

    ngbs_filter = np.zeros(sim.pins.shape[0])
    # Filter rows by checking if any of the elements is less than zero
    # This means that if any of the neighbours of the i-th lattice site is
    # -1, we pin the spin direction at that site
    ngbs_filter = np.any(sim.mesh.neighbours < 0, axis=1, out=ngbs_filter)

    sim.set_pins(ngbs_filter)

# -----------------------------------------------------------------------------


# Debug Information -----------------------------------------------------------

print('-' * 80)
print('Saturation Magnetisation: {} mu_B'.format(args.mu_s))
print('Exchange constant: {}  meV'.format(args.J))
print('DMI constant: {}  meV'.format(args.D))
if args.k_u:
    print('Anisotropy constant: {}   meV'.format(args.k_u))
if args.B:
    print('Zeeman field: (0, 0, {})  T'.format(args.B))
print('-' * 80)

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Initiate simulation ---------------------------------------------------------
# -----------------------------------------------------------------------------


# ================= SIMULATION FUNCTION ======================================

# Define a function to vary the number of initial images and spring constants
"""
Execute a simulation with the NEB function of the fidimag code, for a
stripe (or track or rectangle).

WARNING: interpolation is only well defined in SPHERICAL coordinates,
         thus it is recommended to start the system with the SNEBM
         and later continue the iteration with the CNEBM

vtks and npys are saved in files starting with the 'args.sim_name'
string (only initial and final states)

The frequency of how vtu or npy files are saved is optional

"""

# SIMULATION  ----------------------------------------------------------------

# Interpolations and initial images ------------------------------------------
def load_state(path):
    if os.path.isdir(path):
        if not path.endswith('/'):
            path += '/'
        _file = glob.glob(path + '*.npy')[0]
    else:
        _file = path

    return np.load(_file)

if args.interpolation:
    # Load the states from every 2 arguments of
    # the --interpolation option and the numbers from every  two arguments
    # starting from the 1st element of the list
    images = [load_state(state) for state in args.interpolation[::2]]
    interpolations = [int(num_interps)
                      for num_interps in args.interpolation[1::2]]

    print('==================== Interpolating!  ====================== \n')

elif args.interpolation_rotation:
    images = [load_state(state) for state in args.interpolation_rotation[::2]]
    interpolations = [int(num_interps)
                      for num_interps in args.interpolation_rotation[1::2]]

    print('==================== Interpolating!  ====================== \n')

elif args.images_files:
    # Load the states from the specified npys folder in --images_files
    # If the folder finihes with _LAST, we search the largest number
    # among the folders which we assume finish with a number instead of LAST
    if args.images_files.endswith('_LAST'):
        folders = glob.glob(args.images_files[:-4] + '*')
        folders = sorted(folders,
                         key=lambda f: int(f[len(args.images_files[:-4]):])
                         )

        args.images_files = folders[-1]

    # We will sort the files using the number in their names
    # assuming that they have the structure:
    # 'image_1.npy', 'image_13.npy', etc.
    # For this, we search until two digits with regex and
    # return the integer
    images = [np.load(os.path.join(args.images_files, _file))
              for _file in
              sorted(os.listdir(args.images_files),
                     key=lambda f: int(re.search('\d+', f).group(0))
                     )
              ]

    print('FILES:')
    for _file in sorted(os.listdir(args.images_files),
                        key=lambda f: int(re.search('\d+', f).group(0))
                        ):
        print(_file)

    # No interpolations in this case
    interpolations = None

else:
    print('Specify an accepted method for the initial states')

# ----------------------------------------------------------------------------

# In Fidimag we only use Cartesian coordinates, but the interpolations are
# performed in Spherical coordinates
interp_method = 'linear'
if args.interpolation_rotation:
    interp_method = 'rotation'

if not args.climbing_images:
    args.climbing_images = []

nebm = NEBM_Geodesic(sim,
                     images,
                     interpolations=interpolations,
                     spring_constant=args.nebm_k,
                     name=args.sim_name,
                     climbing_images=args.climbing_images,
                     openmp=True,
                     interpolation_method=interp_method
                     )
if args.tols:
    nebm.create_integrator()
    nebm.set_tols(rtol=args.tols[0], atol=args.tols[1])

# Relax the system with the NEB mehtod
nebm.relax(max_iterations=int(args.nebm_steps),
           save_vtks_every=args.save_vtk,
           save_npys_every=args.save_npy,
           stopping_dYdt=args.stopping_dYdt
           )

if args.save_interpolation:
    if not os.path.exists('dats'):
        os.mkdir('dats')

    # Generate a DAT file with the data from a cubic interpolation for the band
    interp_data = np.zeros((args.save_interpolation, 2))
    (interp_data[:, 0],
     interp_data[:, 1]) = nebm.compute_polynomial_approximation(args.save_interpolation)

    fname = args.sim_name + '_interpolation.dat'
    np.savetxt('dats/' + fname, interp_data)
