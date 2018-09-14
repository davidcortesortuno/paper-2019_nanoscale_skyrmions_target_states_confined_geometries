from __future__ import print_function

"""
"""
# -----------------------------------------------------------------------------

import matplotlib

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np

from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import UniformExchange, DMI, Anisotropy, DemagHexagonal, DemagFull
from fidimag.atomistic import Zeeman
# Import physical constants from fidimag
import fidimag.common.constant as const


import os
import shutil
import types
import re
import sys

this_path = os.path.abspath(__file__)
dir_this_path = os.path.dirname(this_path)

sys.path.append(dir_this_path + '/../../')
import sim_from_image as sfi

sys.path.append(dir_this_path + '/../../mesh_geometries')
from sim_geometries import TruncatedTriangleSim, HexagonSim

# sys.path.append('../../')
# import sim_from_image as sfi


# -----------------------------------------------------------------------------

def hysteresis_loop(config_file,
                    D=1.56, J=5.88, k_u=0.41, mu_s=3, B=(0, 0, 0), Demag=None,
                    ):
    """
    The config file must have the following parameters:
        D
        J
        k_u
        mu_s                :: Magnitude in Bohr magneton units. A file path
                               can be specified to load a NPY file with the
                               mu_s values, when using the mesh_from_image
                               option
        Demag               :: Set to True for Demag
        sim_name            :: Simulation name
        initial_state       :: A function or a npy file
        hysteresis_steps    ::
        mesh_from_image     :: [IMAGE_PATH, xmin, xmax, ymin, ymax]
        hexagonal_mesh      :: [nx, ny, a]
        PBC_2D              :: Set to True for Periodic boundaries
        pin_boundaries      :: Set to True to pin the spins at the boundaries.
                               Their orientations are already given from the
                               initial_state NPY file
        llg_dt              ::
        llg_stopping_dmdt   ::
        llg_max_steps       ::
        llg_do_precession   :: False as default
        llg_alpha           :: 0.01 as default

    """

    # Parameters --------------------------------------------------------------

    cf = {}
    # execfile(config_file, cf)
    # Python3:
    exec(open(config_file).read(), cf)

    D = cf["D"] * const.meV
    J = cf["J"] * const.meV
    k_u = cf["k_u"] * const.meV

    if isinstance(cf["mu_s"], int) or isinstance(cf["mu_s"], float):
        mu_s = cf["mu_s"] * const.mu_B

    if isinstance(cf["initial_state"], str):
        init_state = np.load(cf["initial_state"])
    elif isinstance(cf["initial_state"], types.FunctionType):
        init_state = cf["initial_state"]

    # Set up default arguments
    default_args = {"mesh_alignment": 'diagonal',
                    "mesh_unit_length": 1e-9,
                    "llg_dt": 1e-11,
                    "llg_stopping_dmdt": 1e-2,
                    "llg_max_steps": 4000,
                    "llg_do_precession": False,
                    "llg_alpha": 0.01
                    }

    for key in default_args.keys():
        if not cf.get(key):
            print(default_args[key])
            cf[key] = default_args[key]

    # Simulation object -------------------------------------------------------

    if cf.get("hexagonal_mesh"):
        if not cf["PBC_2D"]:
            mesh = HexagonalMesh(cf["hexagonal_mesh"][2] * 0.5,
                                 int(cf["hexagonal_mesh"][0]),
                                 int(cf["hexagonal_mesh"][1]),
                                 alignment=cf["mesh_alignment"],
                                 unit_length=cf["mesh_unit_length"]
                                 )

        else:
            mesh = HexagonalMesh(cf["hexagonal_mesh"][2] * 0.5,
                                 int(cf["hexagonal_mesh"][0]),
                                 int(cf["hexagonal_mesh"][1]),
                                 periodicity=(True, True),
                                 alignment=cf["mesh_alignment"],
                                 unit_length=cf["mesh_unit_length"]
                                 )

        sim = Sim(mesh, name=cf["sim_name"])

    elif cf.get("mesh_from_image"):
        sim_from_image = sfi.sim_from_image(
            cf["mesh_from_image"][0],
            image_range=[float(cf["mesh_from_image"][1]),
                         float(cf["mesh_from_image"][2]),
                         float(cf["mesh_from_image"][3]),
                         float(cf["mesh_from_image"][4])
                         ],
            sim_name=cf["sim_name"]
            )

        if isinstance(cf["mu_s"], str):
            sim_from_image.generate_magnetic_moments(load_file=cf["mu_s"])
        else:
            sim_from_image.generate_magnetic_moments(mu_s=(mu_s))

        sim = sim_from_image.sim

    elif cf.get("truncated_triangle"):
        if len(cf["truncated_triangle"]) == 3:
            sim_triangle = TruncatedTriangleSim(
                cf["truncated_triangle"][0],  # L
                cf["truncated_triangle"][1],  # offset
                cf["truncated_triangle"][2],  # a
                cf["mu_s"],                   # mu_s
                name=cf["sim_name"]
                )
        elif len(cf["truncated_triangle"]) == 5:
            sim_triangle = TruncatedTriangleSim(
                cf["truncated_triangle"][0],    # L
                [float(offs) for offs in cf["truncated_triangle"][1:4]],  # offsets
                cf["truncated_triangle"][4],    # a
                cf["mu_s"],                     # mu_s
                name=cf["sim_name"]
                )

        sim = sim_triangle.sim

    elif cf.get("hexagon"):
        sim_hexagon = HexagonSim(cf["hexagon"][0],    # R
                                 cf["hexagon"][1],    # a
                                 cf["mu_s"],          # mu_s
                                 name=cf["sim_name"]
                                 )
        sim = sim_hexagon.sim

    # Initial state
    sim.set_m(init_state)

    sim.driver.do_precession = cf["llg_do_precession"]
    sim.driver.alpha = cf["llg_alpha"]

    # Material parameters -----------------------------------------------------

    if cf.get("hexagonal_mesh"):
        sim.mu_s = mu_s

    exch = UniformExchange(J)
    sim.add(exch)

    dmi = DMI(D=(D), dmi_type='interfacial')
    sim.add(dmi)

    zeeman = Zeeman((0, 0, 0))
    sim.add(zeeman, save_field=True)

    if k_u:
        # Uniaxial anisotropy along + z-axis
        sim.add(Anisotropy(k_u, axis=[0, 0, 1]))

    if cf.get("Demag"):
        print('Using Demag!')
        sim.add(DemagHexagonal())

    # Pin boundaries ----------------------------------------------------------

    # We will correct the spin directions according to the specified argument,
    # in case the spins at the boundaries are pinned
    if cf.get('pin_boundaries'):
        ngbs_filter = np.zeros(sim.pins.shape[0])
        # Filter rows by checking if any of the elements is less than zero
        # This means that if any of the neighbours of the i-th lattice site is
        # -1, we pin the spin direction at that site
        ngbs_filter = np.any(sim.mesh.neighbours < 0, axis=1, out=ngbs_filter)

        sim.set_pins(ngbs_filter)

    # Hysteresis --------------------------------------------------------------

    for ext in ['npys', 'vtks']:
        if not os.path.exists('{}/{}'.format(ext, cf["sim_name"])):
            os.makedirs('{}/{}'.format(ext, cf["sim_name"]))

    for ext in ['txts', 'dats']:
        if not os.path.exists('{}/'.format(ext)):
            os.makedirs('{}'.format(ext))

    Brange = cf["hysteresis_steps"]
    print('Computing for Fields:', Brange)

    # We will save the hysteresis steps on this file with every row as:
    # step_number field_in_Tesla
    hystfile = '{}_hyst_steps.dat'.format(cf["sim_name"])
    # If the file already exists, we will append the new steps, otherwise
    # we just create a new file (useful for restarting simulations)
    if not os.path.exists(hystfile):
        nsteps = 0
        fstate = 'w'
    else:
        # Move old txt file from the previous simulation, appending an _r
        # everytime a simulation with the same name is started
        txtfile = [f for f in os.listdir('.') if f.endswith('txt')][0]
        txtfile = re.search(r'.*(?=\.txt)', txtfile).group(0)
        shutil.move(txtfile + '.txt', txtfile + '_r.txt')

        nsteps = len(np.loadtxt(hystfile))
        fstate = 'a'

    f = open(hystfile, fstate)

    for i, B in enumerate(Brange):
        sim.get_interaction('Zeeman').update_field(B)

        sim.driver.relax(dt=cf["llg_dt"],
                         stopping_dmdt=cf["llg_stopping_dmdt"],
                         max_steps=cf["llg_max_steps"],
                         save_m_steps=None,
                         save_vtk_steps=None
                         )

        print('Saving NPY for B = {}'.format(B))
        np.save('npys/{0}/step_{1}.npy'.format(cf["sim_name"], i + nsteps),
                sim.spin)

        sim.driver.save_vtk()
        shutil.move('{}_vtks/m_{}.vtk'.format(cf["sim_name"],
                                              str(sim.driver.step).zfill(6)
                                              ),
                    'vtks/{0}/step_{1}.vtk'.format(cf["sim_name"], i + nsteps)
                    )

        f.write('{} {} {} {}\n'.format(i + nsteps,
                                       B[0], B[1], B[2],
                                       )
                )
        f.flush()
        sim.driver.reset_integrator()

    os.rmdir('{}_vtks'.format(cf["sim_name"]))
    shutil.move('{}.txt'.format(cf["sim_name"]), 'txts/')
    shutil.move(hystfile, 'dats/')

    f.close()

if __name__ == "__main__":

    config_file = sys.argv[1]
    hysteresis_loop(config_file)
