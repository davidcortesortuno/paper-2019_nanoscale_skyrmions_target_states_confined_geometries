from __future__ import print_function
from __future__ import division

"""

Script to run Fidimag simulations using Bash or IPython notebooks
with an optional Real Time View of the dynamics!

These simulations are based on hexagonal lattices.

Since we use argparse, a series of arguments must be specified.

If run from an IPython notebook, use the %run magic, since it
plays nicely with NBagg. For some reason, %%bash does not works with
NBAgg

This program requires Matplotlib >= 1.4.3 to use the Live Preview
in a notebook.

Author: David I. C.
email: d.i.cortes@soton.ac.uk

Modification date: Mon 02 Nov 2015 18:28:03 GMT

"""


# See if this library is called from an IPython notebook and
# return True or False accordingly
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

import sys
import argparse


# -----------------------------------------------------------------------------
# ARGUMENTS -------------------------------------------------------------------
# -----------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='NEB method for 2D films with '
                                 'interfacial DMI')

initial_state = parser.add_mutually_exclusive_group(required=True)

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

parser.add_argument('--D', help='DMI constant in units of meV',
                    type=float, default=1)

parser.add_argument('--J', help='Exchange constant in units of meV',
                    type=float, nargs='+')

parser.add_argument('--mu_s', help='Magnetisation in units of mu_B',
                    default=2)

parser.add_argument('--k_u', help='Anisotropy constant in units of meV',
                    type=float)

parser.add_argument('--B', help='External magnetic field perpendicular to the'
                    ' square plane (z direction), in Tesla',
                    type=float)

parser.add_argument('--Demag', help='Add this option to use dipolar '
                    'interactions using the FFT technique (TESTING)',
                    action='store_true')

parser.add_argument('--FullDemag', help='Add this option to use dipolar '
                    'interactions using a direct calculation (TESTING)',
                    action='store_true')

parser.add_argument('sim_name',
                    help='Simulation name')

parser.add_argument('--PBC_2D',
                    help='Two dimensional boundary condition',
                    action='store_true')

initial_state.add_argument('--initial_state_skyrmion_down',
                           help='This option puts a skyrmionic texture'
                           ' in the centre of the'
                           ' nanotrack, as initial m configuration. The'
                           ' other spins are in the (0, 0, 1) direction',
                           type=float, nargs='+',
                           metavar=('SK_INITIAL_RADIUS')
                           # action='store_true'
                           )

initial_state.add_argument('--initial_state_skyrmion_up',
                           help='This option puts a skyrmionic texture'
                           ' with its core pointing in the +z direction, '
                           'in the centre of the'
                           ' nanotrack, as initial m configuration. The'
                           ' other spins are in the (0, 0, 1) direction',
                           type=float, nargs='+',
                           metavar=('SK_INITIAL_RADIUS')
                           # action='store_true'
                           )

initial_state.add_argument('--initial_state_ferromagnetic_up',
                           help='This option sets the initial '
                           'm configuration as a ferromagnetic state'
                           ' in the (0, 0, 1) direction',
                           action='store_true'
                           )

initial_state.add_argument('--initial_state_ferromagnetic_down',
                           help='This option sets the initial '
                           'm configuration as a ferromagnetic state'
                           ' in the (0, 0, -1) direction',
                           action='store_true'
                           )

initial_state.add_argument('--initial_state_irregular',
                           help='This option sets the initial '
                           'm configuration as an irregular state'
                           ' (TESTING)',
                           action='store_true')

initial_state.add_argument('--initial_state_sk_up_helicoid',
                           help='This option sets the initial '
                           'm configuration as a helicoid with '
                           'a given period (in nm), to '
                           'the right side of the '
                           'mesh, and a skyrmion at a proportion of '
                           'the mesh given '
                           'by the two arguments to this option, '
                           'for the x and y directions (TESTING)',
                           type=float, nargs=3,
                           metavar=('SK_POS_X', 'SK_POS_Y', 'H_PER')
                           )

initial_state.add_argument('--initial_state_sk_down_helicoid',
                           help='',
                           type=float, nargs=3,
                           metavar=('SK_POS_X', 'SK_POS_Y', 'H_PER')
                           )

initial_state.add_argument('--initial_state_sk_up_2helicoid',
                           help='This option sets the initial '
                           'm configuration as a helicoid with '
                           'a given period (in nm), to '
                           'the right side of the '
                           'mesh, and a skyrmion at a proportion of '
                           'the mesh given '
                           'by the two arguments to this option, '
                           'for the x and y directions (TESTING)',
                           type=float, nargs=3,
                           metavar=('SK_POS_X', 'SK_POS_Y', 'H_PER')
                           )

initial_state.add_argument('--initial_state_n_sk_down_helicoid',
                           help='This option sets the initial '
                           'm configuration as a helicoid with '
                           'a given period (in nm), to '
                           'the right side of the '
                           'mesh, and N skyrmions at a proportion of '
                           'the mesh given '
                           'by the two arguments to this option, '
                           'for the x and y directions (TESTING)',
                           type=float, nargs='+',
                           metavar=('SK1_POS_X', 'SK1_POS_Y',
                                    'SKi_POS_X', 'SKi_POS_Y',
                                    'SKN_POS_X', 'SKN_POS_Y',
                                    'H_PER')
                           )

initial_state.add_argument('--initial_state_n_dots_down_helicoid_REL',
                           help='This option sets the initial '
                           'm configuration as a pseudo-helicoid '
                           'at the right side of the '
                           'mesh, and N dots at a proportion of '
                           'the mesh ([-1, 1]) given '
                           'by the two arguments to this option, '
                           'for the x and y directions (TESTING)',
                           type=float, nargs='+',
                           metavar=('DOT1_POS_X', 'DOT1_POS_Y',
                                    'DOTi_POS_X', 'DOTi_POS_Y',
                                    'DOTN_POS_X', 'DOTN_POS_Y',
                                    'HEL_SIZE')
                           )

initial_state.add_argument('--initial_state_multiple_sks_up',
                           help='This option generates multiple skyrmion '
                           'configurations on the sample. The arguments '
                           'are: [sk_radius x_0 y_0 x_1 y_1 ...] '
                           'where x_i,y_i are the coordinates of the '
                           'i-th skyrmion to be generated, given in nm',
                           nargs='+',
                           metavar=('SK_RADIUS', 'SK_POSITIONS')
                           )
initial_state.add_argument('--initial_state_multiple_sks_down',
                           help='This option generates multiple skyrmion '
                           'configurations on the sample. The arguments '
                           'are: [sk_radius x_0 y_0 x_1 y_1 ...] '
                           'where x_i,y_i are the coordinates of the '
                           'i-th skyrmion to be generated, given in nm',
                           nargs='+',
                           metavar=('SK_RADIUS', 'SK_POSITIONS')
                           )

initial_state.add_argument('--initial_state_multiple_dots_down_REL',
                           help='This option generates multiple dots '
                           'in the -z direction. The arguments '
                           'are: [dot_radius(nm) x_0 y_0 x_1 y_1 ...] '
                           'where x_i,y_i are the coordinates of the '
                           'i-th dot, given in relative units in [-1, 1]',
                           nargs='+',
                           metavar=('DOT_RADIUS', 'DOT_POSITIONS')
                           )

initial_state.add_argument('--initial_state_helicoid',
                           help='This option generates a helicoid. The '
                           'options are: helicoid_period helicoid_direction. '
                           'Where the direction is x or y and period is given '
                           'in nm',
                           nargs=2,
                           metavar=('PERIOD', 'DIRECTION')
                           )

initial_state.add_argument('--initial_state_helicoid_rotated',
                           help='This option generates a helicoid. The '
                           'options are: helicoid_period rotation_angle '
                           'Where the period is given in nm and the angle '
                           'in degrees',
                           nargs='+',
                           metavar=('PERIOD', 'DIRECTION_ANGLE')
                           )

initial_state.add_argument('--initial_state_random',
                           help='Randomly oriented spins as initial state',
                           type=int, nargs=1
                           )

initial_state.add_argument('--initial_state_from_image',
                           help='Generate the initial state from the image '
                           'assuming white bg and red/black for the spin '
                           'up/down orientation respectively',
                           action='store_true'
                           )

initial_state.add_argument('--initial_state_file',
                           help='Path to an NPY file with a magnetisation '
                           'field profile'
                           )

parser.add_argument('--preview', help='Specify if instead of relaxing the '
                    'system, it will be shown a real time plot of the '
                    'magnetisation dynamics on the TOP layer (in the z '
                    'direction). This will run for 4 nanoseconds',
                    action='store_true'
                    )

parser.add_argument('--unit_length', help='Mesh unit length',
                    type=float, default=1e-9)

parser.add_argument('--alpha', help='Damping constant value',
                    type=float, default=0.01)

parser.add_argument('--save_files', help='Save vtk and npy files every x'
                    ' steps',
                    type=float, default=None)

parser.add_argument('--save_initial_vtk', help='Save initial VTK step',
                    action='store_true')

parser.add_argument('--stopping_dmdt', help='Specify an specific dm/dt '
                    'threshold value when relaxing the system (default '
                    ' is 0.01)',
                    type=float, default=0.01)

parser.add_argument('--tols', help='Tolerances for the integrator as '
                    'a pair: rtol atol',
                    type=float, nargs=2)

parser.add_argument('--max_steps', help='Specify maximum number of '
                    'steps for the relaxation (default is 5000)',
                    type=int, default=5000)

parser.add_argument('--no_precession', help='To remove LLG precesion term',
                    action='store_true')

parser.add_argument('--alignment', help='Hexagonal lattice alignment '
                    '(square or diagonal)',
                    default='diagonal')

parser.add_argument('--pin_boundaries', help='Set this option to fix the'
                    'spin directions at the boundaries. It must be specified '
                    'the three directions of the spins at this region',
                    nargs=3, type=float)

parser.add_argument('--driver', help='Driver for the minimisation. '
                    ' Default is LLG',
                    default='llg')

parser.add_argument('--save_energy', help='',
                    action='store_true'
                    )

parser.add_argument('--save_Q', help='',
                    action='store_true'
                    )

# Parser arguments
args = parser.parse_args()

# -----------------------------------------------------------------------------

import matplotlib

if args.preview:
    if run_from_ipython():
        matplotlib.use('nbagg')
        print('Using Backend: NBAgg')
    else:
        matplotlib.use('TkAgg')
        print('Using Backend: TkAgg')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np

from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import Exchange, DMI, Anisotropy, DemagHexagonal, DemagFull
from fidimag.atomistic import Zeeman
# Import physical constants from fidimag
import fidimag.common.constant as const

import os
import shutil

this_path = os.path.abspath(__file__)
dir_this_path = os.path.dirname(this_path)

sys.path.append(dir_this_path + '/../../')
import sim_from_image as sfi

sys.path.append(dir_this_path + '/../../mesh_geometries')
from sim_geometries import TruncatedTriangleSim, HexagonSim

mu0 = 4 * np.pi * 1e-7


# -------------------------------------------------------------------------
# Initial states ----------------------------------------------------------
# -------------------------------------------------------------------------

def generate_skyrmion(pos,
                      sim,
                      sign,           # chirality
                      sk_pos=None,    # position
                      sk_radius=1,
                      pi_factor=1.,
                      core_sign=1.,   # -1 for skyrmion down
                      fm_bg_sign=-1.  # ferromagnetic background
                      ):
    """
    Sign will affect the chirality of the skyrmion
    """
    # We will generate a skyrmion in the middle of the stripe
    # with the core pointing down
    if sk_pos is None:
        # Get the x coordinates from the first and last point in the middle
        # of the stripe along the y direction (remember that the
        # mesh indexing is along the x direction)
        mask = sim.mu_s / const.mu_B > 1e-5
        Lmin_x = np.min(sim.mesh.coordinates[:, 0][mask])
        Lmax_x = np.max(sim.mesh.coordinates[:, 0][mask])

        Lmin_y = np.min(sim.mesh.coordinates[:, 1][mask])
        Lmax_y = np.max(sim.mesh.coordinates[:, 1][mask])

        sk_pos = ((Lmax_x - Lmin_x) * 0.5 + Lmin_x,
                  (Lmax_y - Lmin_y) * 0.5 + Lmin_y)

    x = (pos[0] - sk_pos[0])
    y = (pos[1] - sk_pos[1])

    if np.sqrt(x ** 2 + y ** 2) <= sk_radius:
        # Polar coordinates:
        r = (x ** 2 + y ** 2) ** 0.5
        phi = np.arctan2(y, x)
        # This determines the profile we want for the
        # skyrmion
        # Single twisting: k = pi / R
        k = pi_factor * np.pi / sk_radius

        # We define here a 'hedgehog' skyrmion pointing down
        return (sign * np.sin(k * r) * np.cos(phi),
                sign * np.sin(k * r) * np.sin(phi),
                core_sign * np.cos(k * r))
    else:
        return (0, 0, fm_bg_sign)


def generate_skyrmion_up(pos, sim, sign,
                         sk_pos=None,
                         sk_radius=1,
                         pi_factor=1
                         ):

    # We will generate a skyrmion in the middle of the stripe
    # with the core pointing  up
    if sk_pos is None:
        # Get the x coordinates from the first and last point in the middle
        # of the stripe along the y direction (remember that the
        # mesh indexing is along the x direction)
        Lmin = sim.mesh.coordinates[:, 0][sim.mesh.nx * int(sim.mesh.ny * 0.5)]
        Lmax = sim.mesh.coordinates[:, 0][sim.mesh.nx * int(sim.mesh.ny * 0.5)
                                          + sim.mesh.nx - 1]
        sk_pos = ((Lmax - Lmin) * 0.5 + Lmin, sim.mesh.Ly * 0.5)

    x = (pos[0] - sk_pos[0])
    y = (pos[1] - sk_pos[1])

    if np.sqrt(x ** 2 + y ** 2) <= sk_radius:
        # Polar coordinates:
        r = (x ** 2 + y ** 2) ** 0.5
        phi = np.arctan2(y, x)
        # This determines the profile we want for the
        # skyrmion
        # Single twisting: k = pi / R
        k = pi_factor * np.pi / sk_radius

        # We define here a 'hedgehog' skyrmion pointing down
        return (sign * np.sin(k * r) * np.cos(phi),
                sign * np.sin(k * r) * np.sin(phi),
                np.cos(k * r))
    else:
        return (0, 0, -1)


# TODO: Compute the True helicoid wave length ?
def helicoid_m_field(pos, _lambda=3., sign=1,
                     direction='x', h_range=[-1000, 1000]
                     ):
    """
    Generates a helicoid along the x or y direction with a default
    period of 3 nm (chirality according to the DMI)

    h :: range in the x or y direction where the helix will be defined

    """
    if direction == 'x':
        r = pos[0]
    elif direction == 'y':
        r = pos[1]

    if r <= h_range[1] and r >= h_range[0]:
        return (-sign * np.sin(np.pi * r / _lambda),
                0,
                -sign * np.cos(np.pi * r / _lambda))
    else:
        return (0, 0, sign)


def helicoid_m_field_rotated(pos, _lambda=3., sign=-1,
                             rotation_angle=0,
                             reference=(0, 0, 0), h_range=[-2000, 2000]
                             ):
    """
    An helix whose direction is *rotation_angle* with respect to the
    x-direction
    rotation_angle      :: in degrees [0, 360]
    """

    # x, y = pos[0], pos[1]
    T = rotation_angle * np.pi / 180
    u = np.array([np.cos(T), np.sin(T), 0])

    r = np.array(pos)
    r0 = np.array(reference)
    distance = np.dot(r - r0, u)

    if distance >= h_range[0] and distance <= h_range[1]:
        helix = np.array([sign * np.sin(np.pi * np.dot(u, r - r0) / _lambda) * np.cos(T),
                          sign * np.sin(np.pi * np.dot(u, r - r0) / _lambda) * np.sin(T),
                          sign * np.cos(np.pi * np.dot(u, r - r0) / _lambda)])
    else:
        return (0, 0, sign)

    return tuple(helix)


def sk_helicoid_m_field(pos,
                        sim,
                        _lambda=2,
                        sk_pos_prop=(0.3, 0.5),
                        sk_radius=4,
                        sign=1,  # chirality
                        core_sign=1
                        ):
    """
    Generates a skyrmion and a helicoid in the same sample, using
    a period of 2 for the helicoid and the skyrmion located
    at 0.3 times the mesh size along x and 0.5 times along y

    This is currently working for a square shaped mesh or the
    mesh from an image --

    sk_pos_prop     :: Tuple with proportions (from 0 to 1) of the mesh
                       longitudes along the x and y direction where the
                       skyrmion is located

    """
    x, y = pos[0], pos[1]

    if x > 0.58 * sim.mesh.Lx and x < 0.76 * sim.mesh.Lx:
        return helicoid_m_field(pos, _lambda=_lambda, sign=sign)
    else:
        return generate_skyrmion(pos, sim, sign=sign,
                                 sk_pos=(sk_pos_prop[0] * sim.mesh.Lx,
                                         sk_pos_prop[1] * sim.mesh.Ly),
                                 sk_radius=sk_radius,
                                 core_sign=core_sign,
                                 fm_bg_sign=-core_sign
                                 )


def sk_2helicoid_m_field(pos, sim,
                         _lambda=2,
                         sk_pos_prop=(0.4, 0.5),
                         sk_radius=1.3
                         ):
    """
    Generates a skyrmion and a helicoid in the same sample, using
    a period of 2 for the helicoid and the skyrmion located
    at 0.3 times the mesh size along x and 0.5 times along y

    This is currently working for a square shaped mesh or the
    mesh from an image --

    sk_pos_prop     :: Tuple with proportions (from 0 to 1) of the mesh
                       longitudes along the x and y direction where the
                       skyrmion is located

    """
    x, y = pos[0], pos[1]

    if x > 0.58 * sim.mesh.Lx and x < 0.76 * sim.mesh.Lx:
        return helicoid_m_field_rotated(pos, _lambda=_lambda, sign=-1,
                                        rotation_angle=0,
                                        reference=(0.4 * sim.mesh.Lx, 0, 0))
    elif x < 0.43 * sim.mesh.Lx and x > 0.23 * sim.mesh.Lx:
        return helicoid_m_field_rotated(pos, _lambda=_lambda, sign=-1,
                                        rotation_angle=0,
                                        reference=(0.23 * sim.mesh.Lx, 0, 0))
    else:
        return generate_skyrmion(pos, sim, sign=1,
                                 sk_pos=(sk_pos_prop[0] * sim.mesh.Lx,
                                         sk_pos_prop[1] * sim.mesh.Ly),
                                 sk_radius=sk_radius,
                                 core_sign=1.
                                 )


def n_sk_helicoid_m_field(sim,
                          sk_pos_prop=[0.25, 0.25, 0.25, 0.6],
                          _lambda=2,
                          ):
    """
    Generates n skyrmions and a helicoid in the same sample, using
    a period of 2 for the helicoid and the skyrmions located at ...

    This is currently working for a square shaped mesh or the
    mesh from an image --

    sk_pos_prop     :: List with proportions (from 0 to 1) of the mesh
                       longitudes along the x and y direction where the
                       skyrmions are located: x0 y0 x1 y1

    """
    # N SKS ARE GENERATED (FM BACKGROUND):
    skyrmions = sk_pos_prop
    skyrmions[:-1:2] = [s * sim.mesh.Lx for s in skyrmions[:-1:2]]
    skyrmions[1::2] = [s * sim.mesh.Ly for s in skyrmions[1::2]]
    generate_multiple_skyrmions(sign=-1,
                                sk_positions=skyrmions,
                                sk_radius=2,
                                core_sign=-1.,  # -1 for sk down
                                )
    # INSERT A HELIX MANUALLY:
    m_copy = np.copy(sim.spin.reshape(-1, 3))
    for i, pos in enumerate(sim.mesh.coordinates):
        if pos[0] > 0.58 * sim.mesh.Lx and pos[0] < 0.76 * sim.mesh.Lx:
            m_copy[i] = helicoid_m_field_rotated(pos, _lambda=_lambda, sign=1,
                                                 rotation_angle=0,
                                                 reference=(0.4 * sim.mesh.Lx, 0, 0))
    sim.set_m(m_copy.reshape(-1,))


def n_dots_helicoid_m_field(sim,
                            sk_pos_prop=[0.25, 0.25, 0.25, 0.6],
                            hsize=0.1,
                            sk_radius=1.5
                            ):
    """
    Generates n dots and a helicoid in the same sample, using
    a period of 2 for the helicoid and the skyrmions located at ...

    This is currently working for a square shaped mesh or the
    mesh from an image --

    sk_pos_prop     :: List with proportions (from 0 to 1) of the mesh
                       longitudes along the x and y direction where the
                       skyrmions are located: x0 y0 x1 y1

    """
    # N DOTS ARE GENERATED (FM BACKGROUND):
    skyrmions = sk_pos_prop
    skyrmions[:-1:2] = [s * sim.mesh.Lx for s in skyrmions[:-1:2]]
    skyrmions[1::2] = [s * sim.mesh.Ly for s in skyrmions[1::2]]
    generate_multiple_dots(sk_positions=skyrmions,
                           sk_radius=sk_radius,
                           core_sign=-1.,  # -1 for dot down
                           )
    # INSERT A HELIX MANUALLY:
    m_copy = np.copy(sim.spin.reshape(-1, 3))
    for i, pos in enumerate(sim.mesh.coordinates):
        if abs(pos[0] - 0.6 * sim.mesh.Lx) < hsize * sim.mesh.Lx:
            m_copy[i] = (0, 0, -1)
    sim.set_m(m_copy.reshape(-1,))


def generate_multiple_skyrmions(sign,
                                sk_positions,
                                sk_radius,
                                core_sign=1.,   # -1 for skyrmion down
                                ):
    """
    sk_positions        :: A list with [x_0 y_0 x_1 y_1 ...]
                           where (x_i, y_i) is the position of the i-th
                           skyrmion
    """

    def skyrmion_up(pos, index, sign, sk_radius, sk_pos, spins_array,
                    core_sign=1.):
        x, y = (pos[0] - sk_pos[0]), (pos[1] - sk_pos[1])

        if np.sqrt(x ** 2 + y ** 2) <= sk_radius:
            r, phi = (x ** 2 + y ** 2) ** 0.5, np.arctan2(y, x)
            k = np.pi / sk_radius

            return (sign * np.sin(k * r) * np.cos(phi),
                    sign * np.sin(k * r) * np.sin(phi),
                    core_sign * np.cos(k * r))
        else:
            return spins_array.reshape(-1, 3)[index]

    sim.set_m((0, 0, -core_sign))
    sk_positions = np.array([float(x) for x in sk_positions])
    sk_positions.shape = (-1, 2)
    m_copy = np.copy(sim.spin.reshape(-1, 3))
    for sk_pos in sk_positions:
        for i, pos in enumerate(sim.mesh.coordinates):
            m_copy[i] = skyrmion_up(pos, i, sign,
                                    sk_radius,
                                    sk_pos,
                                    m_copy.reshape(-1,),
                                    core_sign=core_sign
                                    )
    sim.set_m(m_copy.reshape(-1,))


def generate_multiple_dots(sk_positions,
                           sk_radius,
                           core_sign=1.,   # -1 for dot down
                           ):
    """
    sk_positions        :: A list with [x_0 y_0 x_1 y_1 ...]
                           where (x_i, y_i) is the position of the i-th
                           skyrmion
    """

    def dot(pos, index, sk_radius, sk_pos, spins_array,
            core_sign=1.):
        x, y = (pos[0] - sk_pos[0]), (pos[1] - sk_pos[1])

        if np.sqrt(x ** 2 + y ** 2) <= sk_radius:

            return (0, 0, core_sign)
        else:
            return spins_array.reshape(-1, 3)[index]

    sim.set_m((0, 0, -core_sign))
    sk_positions = np.array([float(x) for x in sk_positions])
    sk_positions.shape = (-1, 2)
    m_copy = np.copy(sim.spin.reshape(-1, 3))
    for sk_pos in sk_positions:
        for i, pos in enumerate(sim.mesh.coordinates):
            m_copy[i] = dot(pos, i,
                            sk_radius,
                            sk_pos,
                            m_copy.reshape(-1,),
                            core_sign=core_sign
                            )
    sim.set_m(m_copy.reshape(-1,))


def irregular_state(pos):
    m = np.copy(sim.spin)
    m = m.reshape(3, -1).T
    # We will generate a skyrmion in the middle of the stripe
    # (the origin is there) with the core pointing down
    x = (pos[0] - args.hex_nx * args.hex_a * 0.5)
    y = (pos[1] - args.hex_ny * args.hex_a * (3 / 4.) * 0.5)

    if x > 0:
        return (0, 0, 1)
    else:
        return (0, 0, -1)

# -----------------------------------------------------------------------------

# Mesh and Simulation ---------------------------------------------------------

# Check if multiple exchange constants were passed and set the number of
# neighbour shells according to the number of constants
if len(args.J) > 1:
    shells = len(args.J)
else:
    shells = 1

if args.hexagonal_mesh:
    if not args.PBC_2D:
        mesh = HexagonalMesh(args.hexagonal_mesh[2] * 0.5,
                             int(args.hexagonal_mesh[0]),
                             int(args.hexagonal_mesh[1]),
                             alignment=args.alignment,
                             unit_length=args.unit_length,
                             shells=shells
                             )
    else:
        mesh = HexagonalMesh(args.hex_a, args.hex_nx, args.hex_ny,
                             periodicity=(True, True),
                             alignment=args.alignment,
                             unit_length=args.unit_length,
                             shells=shells
                             )
    sim = Sim(mesh, name=args.sim_name, driver=args.driver)

elif args.mesh_from_image:
    sim_from_image = sfi.sim_from_image(args.mesh_from_image[0],
                                        image_range=[float(args.mesh_from_image[1]),
                                                     float(args.mesh_from_image[2]),
                                                     float(args.mesh_from_image[3]),
                                                     float(args.mesh_from_image[4])
                                                     ],
                                        sim_name=args.sim_name,
                                        shells=shells,
                                        driver=args.driver
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
            name=args.sim_name,
            shells=shells,
            driver=args.driver
            )
    elif len(args.truncated_triangle) == 5:
        sim_triangle = TruncatedTriangleSim(
            args.truncated_triangle[0],    # L
            [float(offs) for offs in args.truncated_triangle[1:4]],  # offsets
            args.truncated_triangle[4],    # a
            float(args.mu_s),              # mu_s
            name=args.sim_name,
            shells=shells,
            driver=args.driver
            )

    sim = sim_triangle.sim


elif args.hexagon:
    sim_hexagon = HexagonSim(args.hexagon[0],   # R
                             args.hexagon[1],   # a
                             float(args.mu_s),  # mu_s
                             name=args.sim_name,
                             shells=shells,
                             driver=args.driver
                             )
    sim = sim_hexagon.sim

# The limits of the mesh we can use later
mask = sim.mu_s / const.mu_B > 1e-5
Lmin_x = np.min(sim.mesh.coordinates[:, 0][mask])
Lmax_x = np.max(sim.mesh.coordinates[:, 0][mask])

Lmin_y = np.min(sim.mesh.coordinates[:, 1][mask])
Lmax_y = np.max(sim.mesh.coordinates[:, 1][mask])

dLx = Lmax_x - Lmin_x
dLy = Lmax_y - Lmin_y

# Fidimag simulation options ------------------------------------------------

# sim.set_tols(rtol=1e-10, atol=1e-14)
sim.driver.alpha = args.alpha
# sim.driver.gamma = 1.76e11

if args.no_precession:
    sim.driver.do_precession = False

# Material parameters -------------------------------------------------------

if args.hexagonal_mesh:
    sim.mu_s = float(args.mu_s) * const.mu_B

if shells == 1:
    exch = Exchange(args.J[0] * const.meV)
else:
    exch = Exchange(np.array(args.J) * const.meV)
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
    print('Using Demag!')
    sim.add(DemagHexagonal())

if args.FullDemag:
    print('Using Full Demag calculation!')
    sim.add(DemagFull())

# -------------------------------------------------------------------------


# -----------------------------------------------------------------------------

# Load magnetisation profile ---------------------------------------------

# Change the skyrmion initial configuration according to the
# chirality of the system (give by the DMI constant)
if args.initial_state_skyrmion_down:
    if len(args.initial_state_skyrmion_down) > 1:
        pi_factor = args.initial_state_skyrmion_down[1]
    else:
        pi_factor = 1

    if len(args.initial_state_skyrmion_down) == 3:
        fm_bg = args.initial_state_skyrmion_down[2]
    else:
        fm_bg = 1

    if args.D > 0:
        chirality = -1
    else:
        chirality = 1

    sim.set_m(lambda x: generate_skyrmion(x, sim, chirality,
                                          sk_radius=args.initial_state_skyrmion_down[0],
                                          pi_factor=pi_factor,
                                          core_sign=-1.,
                                          fm_bg_sign=fm_bg
                                          ))

if args.initial_state_skyrmion_up:
    if len(args.initial_state_skyrmion_up) > 1:
        pi_factor = args.initial_state_skyrmion_up[1]
    else:
        pi_factor = 1

    if len(args.initial_state_skyrmion_up) == 3:
        fm_bg = args.initial_state_skyrmion_up[2]
    else:
        fm_bg = -1

    if args.D > 0:
        chirality = 1
    else:
        chirality = -1

    sim.set_m(lambda x: generate_skyrmion(x, sim, chirality,
                                          sk_radius=args.initial_state_skyrmion_up[0],
                                          pi_factor=pi_factor,
                                          core_sign=1.,
                                          fm_bg_sign=fm_bg
                                          )
              )

# The uniform states will be slightly deviated from the z direction
if args.initial_state_ferromagnetic_up:
    sim.set_m((0, 0.9, 0.9))
if args.initial_state_ferromagnetic_down:
    sim.set_m((0, 0.9, -0.9))
if args.initial_state_irregular:
    sim.set_m(irregular_state)
if args.initial_state_helicoid:
    sim.set_m(
        lambda pos: helicoid_m_field(pos,
                                     _lambda=float(args.initial_state_helicoid[0]),
                                     direction=args.initial_state_helicoid[1]
                                     )
        )

if args.initial_state_helicoid_rotated:
    if len(args.initial_state_helicoid_rotated) == 2:
        sim.set_m(
            lambda pos: helicoid_m_field_rotated(pos,
                                                 _lambda=float(args.initial_state_helicoid_rotated[0]),
                                                 rotation_angle=float(args.initial_state_helicoid_rotated[1])
                                                 )
            )
    elif len(args.initial_state_helicoid_rotated) == 5:
        sim.set_m(
            lambda pos: helicoid_m_field_rotated(pos,
                                                 _lambda=float(args.initial_state_helicoid_rotated[0]),
                                                 rotation_angle=float(args.initial_state_helicoid_rotated[1]),
                                                 reference=[float(i) for i in
                                                            args.initial_state_helicoid_rotated[2].split()],
                                                 h_range=[float(args.initial_state_helicoid_rotated[3]),
                                                          float(args.initial_state_helicoid_rotated[4])]
                                                 )
            )

    else:
        raise

if args.initial_state_sk_up_helicoid:
    sim.set_m(
        lambda x: sk_helicoid_m_field(
            x, sim,
            sk_pos_prop=(args.initial_state_sk_up_helicoid[0],
                         args.initial_state_sk_up_helicoid[1]),
            _lambda=args.initial_state_sk_up_helicoid[2]
            )
        )
# Currently working only with D > 0
if args.initial_state_sk_down_helicoid:
    sim.set_m(
        lambda x: sk_helicoid_m_field(
            x, sim,
            sk_pos_prop=(args.initial_state_sk_down_helicoid[0],
                         args.initial_state_sk_down_helicoid[1]),
            _lambda=args.initial_state_sk_down_helicoid[2],
            sign=-1,
            core_sign=-1
            )
        )
if args.initial_state_sk_up_2helicoid:
    sim.set_m(
        lambda x: sk_2helicoid_m_field(
            x, sim,
            sk_pos_prop=(args.initial_state_sk_up_2helicoid[0],
                         args.initial_state_sk_up_2helicoid[1]),
            _lambda=args.initial_state_sk_up_2helicoid[2]
            )
        )
if args.initial_state_n_sk_down_helicoid:
    n_sk_helicoid_m_field(sim,
                          sk_pos_prop=[skpos for skpos
                                       in args.initial_state_n_sk_down_helicoid[:-1]
                                       ],
                          # Last argument is helix period
                          _lambda=args.initial_state_n_sk_down_helicoid[-1]
                          )
if args.initial_state_n_dots_down_helicoid_REL:
    rad = args.initial_state_n_dots_down_helicoid_REL[-2] * dLx
    n_dots_helicoid_m_field(sim,
                            sk_pos_prop=[skpos for skpos
                                         in args.initial_state_n_dots_down_helicoid_REL[:-2]
                                         ],
                            # Last argument is helix period
                            hsize=args.initial_state_n_dots_down_helicoid_REL[-1],
                            sk_radius=rad
                            )
if args.initial_state_multiple_sks_up:
    if args.D > 0:
        chirality = 1
    else:
        chirality = -1

    generate_multiple_skyrmions(
        sign=chirality,
        sk_positions=args.initial_state_multiple_sks_up[1:],
        sk_radius=float(args.initial_state_multiple_sks_up[0]),
        core_sign=1.
        )
if args.initial_state_multiple_sks_down:
    if args.D > 0:
        chirality = -1
    else:
        chirality = 1
    generate_multiple_skyrmions(
        sign=chirality,
        sk_positions=args.initial_state_multiple_sks_down[1:],
        sk_radius=float(args.initial_state_multiple_sks_down[0]),
        core_sign=-1
        )

if args.initial_state_multiple_dots_down_REL:
    positions = np.array([float(x) for x in args.initial_state_multiple_dots_down_REL[1:]])
    positions[::2] = (positions[::2] + 1) * 0.5 * dLx
    positions[1::2] = (positions[1::2] + 1) * 0.5 * dLy
    print(positions)

    # Radius as a proportion of lenght in the x direction
    rad = float(args.initial_state_multiple_dots_down_REL[0]) * dLx

    generate_multiple_dots(
        sk_positions=positions,
        sk_radius=rad,
        core_sign=-1
        )

if args.initial_state_file:
    sim.set_m(np.load(args.initial_state_file))

if args.initial_state_random:
    np.random.seed = args.initial_state_random
    sim.set_m(lambda r: np.random.uniform(-1, 1, 3))

if args.initial_state_from_image:
    sim_from_image.generate_magnetic_configuration(
        orientations=[(0, 0, -1), (0, 0, 1)],
        colours=[(0., 0., 0.), (1., 0., 0.)]
        )

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
# -------------------------------------------------------------------------


# Debug Information -------------------------------------------------------
# print 'Simulating a {} x {} x {} box'.format(args.box_length,
#                                              args.box_width,
#                                              args.box_thickness)

print('Saturation Magnetisation: {} mu_B'.format(args.mu_s))
print('Exchange constant: {}  meV'.format(args.J))
print('DMI constant: {}  meV'.format(args.D))
if args.k_u:
    print('Anisotropy constant: {}   meV'.format(args.k_u))
if args.B:
    print('Zeeman field: (0, 0, {})  T'.format(args.B))
print('------------------------------------------------------------------')

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Initiate simulation -----------------------------------------------------
# -------------------------------------------------------------------------

# Here we use matplotlib in interactive mode, updating the figure by
# drawing the canvas instead of creating a whole new plot
if args.preview:

    times = np.linspace(0, 4e-9, 500)

    # Data for the Figure
    # Mesh coordinates
    coords = sim.mesh.coordinates
    # Now we copy the magnetisation to not modify the simulation object
    m = np.copy(sim.spin)
    # Reshape and transpose to get an array
    # with [mx, my, mz] elements (matrix)
    m = m.reshape(-1, 3)

    # Figure specifications -----------------------------------------------

    # Aspect according to the square dimensions
    # aspect = args.hex_nx / float(args.hex_ny)
    aspect = 1
    w, h = plt.figaspect(aspect)
    fig = plt.figure(figsize=(w, h))

    ax = fig.add_subplot(111)

    # Colourise the spins according to the mz component
    quiv = ax.quiver(coords[:, 0], coords[:, 1],
                     m[:, 0], m[:, 1], m[:, 2],
                     # Arrow properties (can vary according to the plot)
                     cmap='RdYlBu', width=.008, linewidth=1,
                     scale=1 / 0.05,
                     # Data limits for the colour map
                     clim=[-1, 1]
                     )

    ttime = ax.text(1., 1.05, '',
                    transform=ax.transAxes,
                    # Vertical and horizontal alignment
                    va='center', ha='right')

    tenergy = ax.text(0, 1.05, '',
                      transform=ax.transAxes,
                      va='center', ha='left')

    # Colour bar (append to not distort the main plot)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)

    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    cbar = matplotlib.colorbar.ColorbarBase(cax,
                                            cmap='RdYlBu',
                                            norm=norm,
                                            ticks=[-1, 0, 1],
                                            orientation='vertical',
                                            )

    cbar.set_label(r'$m_z$', rotation=270, labelpad=10, fontsize=16)

    # Interactive mode (this needs so set up a proper backend
    # when importing matplotlib for the first time)
    plt.ion()
    # Set False to avoid the execution of the following code
    plt.show(False)

    # ---------------------------------------------------------------------

    sim.save_vtk()
    # Now run the simulation printing the energy
    for time in times:
        if not run_from_ipython():
            print('Time: ', time, ' s')
            print('Total energy: ', sim.compute_energy(), ' J')
            print('\n')
        sim.driver.run_until(time)

        # Update the vector data for the plot (the spins do not move
        # so we don't need to update the coordinates) and redraw
        m = np.copy(sim.spin)
        # reshape rows, transpose and filter according to top layer
        m = m.reshape(-1, 3)
        quiv.set_UVC(m[:, 0], m[:, 1], m[:, 2])

        # Update title
        ttime.set_text('Time: {:.4f} ns'.format(time * 1e9))
        tenergy.set_text('Energy: {:.6e} J'.format(sim.compute_energy()))

        # fig.show()
        fig.canvas.draw()

    sim.save_vtk()

else:
    # print(np.min(sim.mesh.coordinates[:, 1]), np.max(sim.mesh.coordinates[:, 1]))
    if args.save_initial_vtk:
        sim.save_vtk()
    # Fidimag automatically saves the last state (?)
    sim.driver.do_precession = False

    if args.tols:
        sim.set_tols(rtol=args.tols[0], atol=args.tols[1])

    if args.driver == 'llg':
        sim.driver.relax(dt=1e-13, stopping_dmdt=args.stopping_dmdt,
                         max_steps=args.max_steps,
                         save_m_steps=args.save_files,
                         save_vtk_steps=args.save_files)

        # Set OpenMP for the integrator
        sim.driver.set_integrator('sundials_openmp', use_jac=False)

    elif args.driver == 'minimiser':
        sim.driver.minimise(max_steps=args.max_steps,
                            stopping_dm=args.stopping_dmdt,
                            save_data_steps=1000)

    elif args.driver == 'steepest_descent':
        sim.driver.tmax = 1e-2
        sim.driver.tmin = 1e-16

        sim.driver.minimise(max_steps=args.max_steps,
                            stopping_dm=args.stopping_dmdt,
                            save_data_steps=1000,
                            save_m_steps=args.save_files,
                            save_vtk_steps=args.save_files
                            )

    # Save final states
    sim.save_m()
    sim.save_vtk()


# -------------------------------------------------------------------------
# Files -------------------------------------------------------------------
# -------------------------------------------------------------------------

npy_dir = 'npys/'
vtk_dir = 'vtks/'
txt_dir = 'txts/'

if not os.path.exists(npy_dir):
    os.makedirs(npy_dir)
if not os.path.exists(vtk_dir):
    os.makedirs(vtk_dir)
if not os.path.exists(txt_dir):
    os.makedirs(txt_dir)

if args.save_energy:
    energies_dir = 'energies/'
    if not os.path.exists(energies_dir):
        os.makedirs(energies_dir)

    E = sim.compute_energy()
    np.savetxt(energies_dir + args.sim_name + '_E.txt', np.array([E]))

if args.save_Q:
    Q_dir = 'top_charges/'
    if not os.path.exists(Q_dir):
        os.makedirs(Q_dir)

    Q = sim.skyrmion_number(method='BergLuscher')
    np.savetxt(Q_dir + args.sim_name + '_Q.txt', np.array([Q]))

# files = [_f for _f in os.listdir('.')
#          if (os.path.isdir(_f) and _f.startswith(args.sim_name))]

# Fidimag vtk and npy files are saved in the sim_name_vtks or sim_name_npys
# folders respectively. We will move them to the vtks/ and npys/ directories

# Remove the vtk folder if it already exists inside the vtks/ folder
if os.path.exists(vtk_dir + args.sim_name + '_vtks/'):
    shutil.rmtree(vtk_dir + args.sim_name + '_vtks/')
try:
    shutil.move(args.sim_name + '_vtks/', vtk_dir)
except IOError:
    pass

# Same for npy folder
if os.path.exists(npy_dir + args.sim_name + '_npys/'):
    shutil.rmtree(npy_dir + args.sim_name + '_npys/')
try:
    shutil.move(args.sim_name + '_npys/', npy_dir)
except IOError:
    pass

# Now do this for the txt files
if os.path.exists(txt_dir + args.sim_name + '.txt'):
    os.remove(txt_dir + args.sim_name + '.txt')
try:
    shutil.move(args.sim_name + '.txt', txt_dir)
except IOError:
    pass


