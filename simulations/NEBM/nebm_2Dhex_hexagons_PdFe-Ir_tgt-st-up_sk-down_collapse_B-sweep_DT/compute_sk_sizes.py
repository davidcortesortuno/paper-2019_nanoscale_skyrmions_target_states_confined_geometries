import sys, os, shutil
import numpy as np
import fidimag.common.constant as C
from fidimag.atomistic import UniformExchange, DMI, Anisotropy
from fidimag.atomistic import Zeeman

this_path = os.path.abspath(__file__)
dir_this_path = os.path.dirname(this_path)
sys.path.append(dir_this_path + '/../../../mesh_geometries')
from sim_geometries import TruncatedTriangleSim, HexagonSim

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
B_range = range(0, 1101, 100)
R_range = range(6, 15, 2)
hexagons_tgt_st_m = {}

for i, R in enumerate(R_range):

    hexagons_tgt_st_m[R] = {}

    for j, B in enumerate(B_range):

        if R == 6 and B >= 600:
            continue

        hexagons_tgt_st_m[R][B] = {'x': [], 'mz': []}

        fname = "2Dhex_hexagon_R{}nm_PdFe-Ir_1-tgt-st-up_B{}mT_npys/".format(R, B)
        npy_file = os.listdir(base_folder + fname)[0]

        sim = create_simulation(R, B * 1e-3)
        sim.set_m(np.load(base_folder + fname + npy_file))

        x, y = sim.mesh.coordinates[:, 0], sim.mesh.coordinates[:, 1]
        xs, ys = np.unique(x), np.unique(y)
        mask = np.logical_and(y == ys[int(ys.shape[0] * 0.5)],
                              sim.mu_s / C.mu_B > 1e-6
                              )

        centre_x = (np.max(x[mask]) - np.min(x[mask])) * 0.5 + np.min(x[mask])

        hexagons_tgt_st_m[R][B]['x'] = x[mask] - centre_x
        hexagons_tgt_st_m[R][B]['mz'] = np.copy(sim.spin.reshape(-1, 3)[:, 2][mask])

# -----------------------------------------------------------------------------

st_sizes_mz = {}
st_r_ring = {}

xoffset = 0.1

for i, R in enumerate(R_range):
    st_sizes_mz[R] = {}
    st_r_ring[R] = {}
    for j, B in enumerate(B_range):

        if R == 6 and B >= 600:
            continue

        xmin, xmax = np.min(hexagons_tgt_st_m[R][B]['x']), np.max(hexagons_tgt_st_m[R][B]['x'])

        x_interp = np.linspace(xmin, xmax, 1000)
        f = si.interp1d(hexagons_tgt_st_m[R][B]['x'],
                        hexagons_tgt_st_m[R][B]['mz'], kind='cubic')
        # y_interp = f(x_interp)

        # Inner target state radius radius: m_z = 0
        # Find the minima
        x_minimum = so.fmin(f, 0 + xoffset, disp=False)

        # We compute m_z = 0, starting from x=0 until x_minimum
        # r_tgt_st = so.brentq(f, 0, x_minimum)
        # st_sizes_mz[R][B] = r_tgt_st
        st_r_ring[R][B] = x_minimum

# print(st_sizes_mz[8])


# Generate initial states -----------------------------------------------------

import scipy.interpolate


def move_state(n, n_images, sim_bg, R, B):
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
    # We will use this interpolation function to extract the spin components
    # at any coordinate. This interpolator works better than a
    # nearest neighbour interpolator
    st_field = []
    for i in range(3):
        st_field.append(scipy.interpolate.CloughTocher2DInterpolator(
            sim_bg.mesh.coordinates[:, :2],
            sim_bg.spin.reshape(-1, 3)[:, i]
            ))

    # Mask where ther eis material, which defines the hecagonal island
    mask = (sim_bg.mu_s / C.mu_B) > 1e-5

    # Get an approximated sk radius from the dictionary defined at the
    # beginning of the script
    outer_ring = st_r_ring[R][B]
    inner_ring = outer_ring - (n * outer_ring / n_images)

    # Copy the spin field from the background simulation. We will use
    # this as a canvas to draw the new state
    interp_state = np.copy(sim_bg.spin.reshape(-1, 3))

    # Compute the coordinates of the centre of the hexagon
    x, y = sim_bg.mesh.coordinates[:, 0], sim_bg.mesh.coordinates[:, 1]
    xs, ys = np.unique(x), np.unique(y)
    mask = np.logical_and(y == ys[int(ys.shape[0] * 0.5)],
                          sim_bg.mu_s / C.mu_B > 1e-6)

    centre_x = (np.max(x[mask]) - np.min(x[mask])) * 0.5 + np.min(x[mask])
    centre_y = (np.max(y[mask]) - np.min(y[mask])) * 0.5 + np.min(y[mask])

    # Loop through the coordinates of the lattice
    for i, pos in enumerate(sim_bg.mesh.coordinates):
        # Only process lattice sites with material
        if (sim_bg.mu_s[i] / C.mu_B) > 1e-5:
            xrel, yrel = pos[0] - centre_x, pos[1] - centre_y

            if (xrel ** 2 + yrel ** 2 <= outer_ring ** 2) and \
                    (xrel ** 2 + yrel ** 2 >= inner_ring ** 2):
                # We use the x, y offsets from the new sk to obtain the
                # spin orientation from the original sk, which is
                # positioned at xc, yc
                interp_state[i] = [0, 0, -1]

    return interp_state.reshape(-1,)


N_IMAGES = 26

for R in R_range:
    for B in B_range:

        if R == 6 and B >= 600:
            continue

        # Create directories ..................................................

        npy_folder = 'initial_states/npys/' + 'R{}nm_B{}mT/'.format(R, B)
        vtk_folder = 'initial_states/vtks/' + 'R{}nm_B{}mT/'.format(R, B)
        for _dir in [npy_folder, vtk_folder]:
            if not os.path.exists(_dir):
                os.makedirs(_dir)

        # Set the simulation with the skyrmion and save as the 0th image
        fname_st = "2Dhex_hexagon_R{}nm_PdFe-Ir_1-sk-down_B{}mT_npys/".format(R, B)
        npy_file_st = os.listdir(base_folder + fname_st)[0]
        sim = create_simulation(R, B * 1e-3)
        sim.set_m(np.load(base_folder + fname_st + npy_file_st))
        np.save(npy_folder + 'image_{:06}.npy'.format(N_IMAGES + 1), sim.spin)
        sim.save_vtk()
        shutil.move('unnamed_vtks/m_000000.vtk',
                    vtk_folder + 'image_{:06}.vtk'.format(N_IMAGES + 1))

        # Set the simulation with the FM state and save as the last image
        fname_bg = "2Dhex_hexagon_R{}nm_PdFe-Ir_1-tgt-st-up_B{}mT_npys/".format(R, B)
        npy_file_bg = os.listdir(base_folder + fname_bg)[0]
        sim_bg = create_simulation(R, B * 1e-3)
        sim_bg.set_m(np.load(base_folder + fname_bg + npy_file_bg))
        np.save(npy_folder + 'image_{:06}.npy'.format(0), sim_bg.spin)
        sim_bg.save_vtk()
        shutil.move('unnamed_vtks/m_000000.vtk',
                    vtk_folder + 'image_{:06}.vtk'.format(0))

        # .....................................................................
        # Here we generate the states with the skyrmion shifted in the +y
        # direction

        for i in range(1, N_IMAGES + 1):

            # Set the simulation with the sk at the centre of the hexagon
            sim.set_m(np.load(base_folder + fname_st + npy_file_st))

            interp_state = move_state(i, N_IMAGES, sim_bg, R, B)
            sim.set_m(interp_state)
            # Save in the corresponding folder, for the R, B simulation
            np.save(npy_folder + 'image_{:06}.npy'.format(i), sim.spin)
            sim.save_vtk()
            shutil.move('unnamed_vtks/m_000000.vtk',
                        vtk_folder + 'image_{:06}.vtk'.format(i))
            shutil.rmtree('unnamed_vtks')
