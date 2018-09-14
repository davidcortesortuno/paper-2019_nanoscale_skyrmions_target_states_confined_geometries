import sys, os, shutil
import numpy as np
import fidimag.common.constant as C
from fidimag.atomistic import UniformExchange, DMI, Anisotropy
from fidimag.atomistic import Zeeman

sys.path.append('../../../mesh_geometries')
from sim_geometries import HexagonSim

import scipy.interpolate as si
import scipy.optimize as so


# -----------------------------------------------------------------------------

def create_simulation(R, B):
    mu_s = 3
    sim_hexagon = HexagonSim(R,       # R
                             0.2715,  # a
                             mu_s,    # mu_s
                             name='unnamed'
                             )
    sim = sim_hexagon.sim

    # mask = (sim.mu_s / C.mu_B) > 1e-5
    exch = UniformExchange(5.881 * C.meV)
    sim.add(exch)
    dmi = DMI(D=1.557 * C.meV, dmi_type='interfacial')
    sim.add(dmi)
    sim.add(Zeeman([0., 0., B]))
    sim.add(Anisotropy(0.406 * C.meV, axis=[0, 0, 1]))

    return sim

# -----------------------------------------------------------------------------

base_folder = '../../relaxation/hexagons_size_variation_DT/npys/'
B_range = range(0, 2001, 100)
hexagons_sk_down_m = {}

for i, R in enumerate([3, 4, 6, 8, 10, 12]):

    hexagons_sk_down_m[R] = {}

    for j, B in enumerate(B_range):

        hexagons_sk_down_m[R][B] = {'x': [], 'mz': []}

        fname = "2Dhex_hexagon_R{}nm_PdFe-Ir_1-sk-down_B{}mT_npys/".format(R, B)
        npy_file = os.listdir(base_folder + fname)[0]

        sim = create_simulation(R, B * 1e-3)
        sim.set_m(np.load(base_folder + fname + npy_file))

        x, y = sim.mesh.coordinates[:, 0], sim.mesh.coordinates[:, 1]
        xs, ys = np.unique(x), np.unique(y)
        mask = np.logical_and(y == ys[int(ys.shape[0] * 0.5)],
                              sim.mu_s / C.mu_B  > 1e-6
                              )

        centre_x = (np.max(x[mask]) - np.min(x[mask])) * 0.5 + np.min(x[mask])

        hexagons_sk_down_m[R][B]['x'] = x[mask] - centre_x
        hexagons_sk_down_m[R][B]['mz'] = np.copy(sim.spin.reshape(-1, 3)[:, 2][mask])

# -----------------------------------------------------------------------------

sk_sizes_mz = {}

for R in [10, 8, 6, 4, 3]:
    sk_sizes_mz[R] = {}
    xmin, xmax = np.min(hexagons_sk_down_m[R][0]['x']), np.max(hexagons_sk_down_m[R][0]['x'])

    if R == 10:
        B_range = range(700, 2001, 100)
    else:
        B_range = range(0, 2001, 100)

    for i, B in enumerate(B_range):

        x_interp = np.linspace(xmin, xmax, 1000)
        f = si.interp1d(hexagons_sk_down_m[R][B]['x'],
                        hexagons_sk_down_m[R][B]['mz'], kind='cubic')
        # y_interp = f(x_interp)

        # Skyrmion radius: m_z = 0
        r_sk = so.brentq(f, 0, xmax)

        sk_sizes_mz[R][B] = r_sk

# print(sk_sizes_mz[8])


# Generate initial states -----------------------------------------------------

import scipy.interpolate


def move_skyrmion(n, n_images, sim_sk, sim_fm, R, B):
    """
    From two simulations: one where a spin is at the middle of the hexagon,
    sim_sk, and one with the uniform state, sim_fm; this function generates a
    spin field array where a skyrmion is shifted towards the +y direction
    by a distance given by *n*.

    The new skyrmion y-position is computed as:

        hexagon_centre + n * (max_y + sk_radius - hexagon_centre) / n_images

    where max_y is the y coordinate at the top of the hexagon. This means
    that we divide the length between the hexagon centre and the top of
    the hexagon in n_images points. Then we position the skyrmion at the
    n-th point. We add the sk_radius quantity so the skyrmion completely
    escapes from the hexagon for n close to n_images (if don't, half of the
    skyrmion stays at the sample at the latest point)

    R, B: radius of the hexagon and magnetic field magnitude
          (to use the right simulation)

    """

    # Interpolate the spin field from the simulation with the skyrmion
    # We will use this interpolation function to create the skyrmion at
    # the shifted position. This interpolator works better than a
    # nearest neighbour interpolator
    sk_field = []
    for i in range(3):
        sk_field.append(scipy.interpolate.CloughTocher2DInterpolator(
            sim_sk.mesh.coordinates[:, :2],
            sim_sk.spin.reshape(-1, 3)[:, i]
            ))

    # Mask where ther eis material, which defines the hecagonal island
    mask = (sim_sk.mu_s / C.mu_B) > 1e-5

    # Compute the coordinates of the centre of the hexagon
    rs = sim_sk.mesh.coordinates
    xs, ys = np.unique(rs[:, 0][mask]), np.unique(rs[:, 1][mask])
    xc, yc = xs[int(len(xs) * 0.5)], ys[int(len(ys) * 0.5)]

    # Get an approximated sk radius from the dictionary defined at the
    # beginning of the script. We add 1 since the approximation m_z = 0
    # is not very accurate for the skyrmion size
    sk_radius = sk_sizes_mz[R][B] + 1

    # The new position for the skyrmion that we will generate
    shifted_sk_pos = (xc, yc + n * (sk_radius + ys.max() - yc) / n_images)
    # Copy the spin field from the uniform state simulation. We will use
    # this as a canvas to draw the new skyrmion with shifted position
    interp_state = np.copy(sim_fm.spin.reshape(-1, 3))

    # Loop through the coordinates of the lattice
    for i, pos in enumerate(sim_fm.mesh.coordinates):
        # Only process lattice sites with material
        if (sim_fm.mu_s[i] / C.mu_B) > 1e-5:
            # Position relative to the centre of the sk with shifted position
            xrel, yrel = pos[0] - shifted_sk_pos[0], pos[1] - shifted_sk_pos[1]
            # If it is smaller than the original sk radius, copy the
            # corresponding spin orientation from the sim with the skyrmion
            # We obtain the spin direction from the interpolation function
            # for the sk, which allows us to pass a coordinate
            if xrel ** 2 + yrel ** 2 <= sk_radius ** 2:
                # We use the x, y offsets from the new sk to obtain the
                # spin orientation from the original sk, which is
                # positioned at xc, yc
                interp_state[i] = [sk_field[0](xc + xrel, yc + yrel),
                                   sk_field[1](xc + xrel, yc + yrel),
                                   sk_field[2](xc + xrel, yc + yrel)]

    return interp_state.reshape(-1,)


N_IMAGES = 26

for R in [10, 8, 6, 4, 3]:

    if R == 10:
        B_range = range(700, 2001, 100)
    else:
        B_range = range(0, 2001, 100)

    for B in B_range:

        # Create directories ..................................................

        npy_folder = 'initial_states/npys/' + 'R{}nm_B{}mT/'.format(R, B)
        vtk_folder = 'initial_states/vtks/' + 'R{}nm_B{}mT/'.format(R, B)
        for _dir in [npy_folder, vtk_folder]:
            if not os.path.exists(_dir):
                os.makedirs(_dir)

        # Set the simulation with the skyrmion and save as the 0th image
        fname_sk = "2Dhex_hexagon_R{}nm_PdFe-Ir_1-sk-down_B{}mT_npys/".format(R, B)
        npy_file_sk = os.listdir(base_folder + fname_sk)[0]
        sim = create_simulation(R, B * 1e-3)
        sim.set_m(np.load(base_folder + fname_sk + npy_file_sk))
        np.save(npy_folder + 'image_{:06}.npy'.format(0), sim.spin)
        sim.save_vtk()
        shutil.move('unnamed_vtks/m_000000.vtk',
                    vtk_folder + 'image_{:06}.vtk'.format(0))

        # Set the simulation with the FM state and save as the last image
        fname_fm = "2Dhex_hexagon_R{}nm_PdFe-Ir_fm-up_B{}mT_npys/".format(R, B)
        npy_file_fm = os.listdir(base_folder + fname_fm)[0]
        sim_fm = create_simulation(R, B * 1e-3)
        sim_fm.set_m(np.load(base_folder + fname_fm + npy_file_fm))
        np.save(npy_folder + 'image_{:06}.npy'.format(N_IMAGES + 1), sim_fm.spin)
        sim_fm.save_vtk()
        shutil.move('unnamed_vtks/m_000000.vtk',
                    vtk_folder + 'image_{:06}.vtk'.format(N_IMAGES + 1))

        # .....................................................................
        # Here we generate the states with the skyrmion shifted in the +y
        # direction

        for i in range(1, N_IMAGES + 1):

            # Set the simulation with the sk at the centre of the hexagon
            sim.set_m(np.load(base_folder + fname_sk + npy_file_sk))

            # Get the spin field from the shifted skyrmion and set this
            # new state in the sk simulation (or the FM simulation, just
            # dont forget to reload it at the beginning)
            interp_state = move_skyrmion(i, N_IMAGES, sim, sim_fm, R, B)
            sim.set_m(interp_state)
            # Save in the corresponding folder, for the R, B simulation
            np.save(npy_folder + 'image_{:06}.npy'.format(i), sim.spin)
            sim.save_vtk()
            shutil.move('unnamed_vtks/m_000000.vtk',
                        vtk_folder + 'image_{:06}.vtk'.format(i))
            shutil.rmtree('unnamed_vtks')


