from __future__ import print_function

# -----------------------------------------------------------------------------
import ipywidgets as ipw
from IPython.display import display

import fidimag
from fidimag.common import CuboidMesh
from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
import fidimag.common.constant as const

from fidimag.atomistic import UniformExchange, DMI, Anisotropy, DemagHexagonal
from fidimag.atomistic import Zeeman

import sim_from_image as sfi

import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from datetime import datetime, timedelta
from functools import wraps

import colorsys

# -----------------------------------------------------------------------------


# Taken from https://gist.github.com/ChrisTM/5834503
class throttle(object):
    """
    Decorator that prevents a function from being called more than once every
    time period.
    To create a function that cannot be called more than once a minute:
        @throttle(minutes=1)
        def my_fun():
            pass
    """
    def __init__(self, seconds=0, minutes=0, hours=0):
        self.throttle_period = timedelta(
            seconds=seconds, minutes=minutes, hours=hours
        )
        self.time_of_last_call = datetime.min

    def __call__(self, fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            now = datetime.now()
            time_since_last_call = now - self.time_of_last_call

            if time_since_last_call > self.throttle_period:
                self.time_of_last_call = now
                return fn(*args, **kwargs)

        return wrapper


# -----------------------------------------------------------------------------


class interactive_simulation(object):

    def __init__(self, config_file=False,
                 simulation='2D_square',
                 D=1.56, J=5.88, ku=0.41, mu_s=3, B=(0, 0, 0), Demag=None,
                 mesh_nx=50, mesh_ny=50, mesh_a=0.2715
                 ):
        """

        This class generates an interactive simulation using Ipython widgets.
        By default, the image is a square shaped hexagonal mesh with finite
        boundaries.

        Arguments:

        config_file        :: Not implemented yet

        simulation         :: Any of the following strings:
                              ['2D_square', 'experimental_sample', '1D_chain']

                              The 2D_square option creates a square shaped
                              mesh using an atomistic hexagonal lattice.

                              The experimental_sample option creates a mesh
                              using an image from experimental data and
                              populates it with spins using the atomistic
                              hexagonal mesh. Here the mesh parameters are not
                              used in the class.  The image processing is done
                              through the sim_from_image library.

                              The 1D_chain creates a one dimensional spin chain

        Default magnetic parameters are from Romming et al. [PRL 114, 177203
        (2015)] research on PdFe bilayers on Ir(111)

        Basic usage:

            # Simulation with default magnetic field along z of 0 T
            i_sim = interactive_simulation()

            i_sim.generate_simulation()
            i_sim.generate_m_field()

            # This produces the plot and simulation options
            i_sim.generate_widgets()


        TODO: Add option to pass a custom function for the magnetisation field
              initial state

              Color plot according to magnetisation directions in 360 degrees?

        """

        self.simulation = simulation

        if config_file:
            tmp_config = {}
            configs = execfile(config_file, tmp_config)

            self.D = configs["D"] * const.meV
            self.J = configs["J"] * const.meV
            self.ku = configs["ku"] * const.meV
            self.mu_s = configs["mu_s"] * const.mu_B
            self.m_field = configs["m_field"]
            if configs["B"] is not None:
                self.B = configs["B"]

        else:
            self.D = D * const.meV
            self.J = J * const.meV
            self.ku = ku * const.meV
            self.mu_s = mu_s * const.mu_B
            self.B = B
            self.Demag = Demag

            self.mesh_nx = mesh_nx
            self.mesh_ny = mesh_ny
            self.mesh_a = mesh_a

        # Dictionary to translate a vector component into the corresponding
        # indexes in Fidimag arrays, i.e. x --> 0, y --> 1, z --> 2
        self.v_dict = {'x': 0, 'y': 1, 'z': 2}

        # Measure for dm / dt
        self.DEGREE_PER_NANOSECOND = 2 * np.pi / (360 * 1e-9)


    def skyrmion_m_field(self, pos, sign,
                         sk_pos=None, sk_r=4, core=1, pi_factor=1.,
                         out_skyrmion_dir=None
                         ):
        """

        Generate a skyrmionic configuration (an approximation) of radius sk_r
        at the (sk_pos[0], sk_pos[1]) position of the mesh

        core      :: Tells if skyrmion is up or down, it can be +1 or -1
        sign      :: Changes the sk chirality, it can be +1 or -1
        pi_factor :: The skyrmion field is generated by the formula:
                          (sign * np.sin(k * r) * np.cos(phi), ...)
                     where k = pi_factor * np_pi / sk_r. The factor determines
                     the number of twistings of the skyrmion profile
        out_skyrmion_dir :: Direction outside the skyrmion (by default it is
                            the direction opposite to the skyrmion core)

        """

        if sk_pos is None:
            # We assume a square sized hexagonal mesh so the centre
            # is at half of every dimension
            sk_pos = self.sim.mesh.Lx * 0.5, self.sim.mesh.Ly * 0.5

        x = (pos[0] - sk_pos[0])
        y = (pos[1] - sk_pos[1])

        if np.sqrt(x ** 2 + y ** 2) <= sk_r:
            # Polar coordinates:
            r = (x ** 2 + y ** 2) ** 0.5
            phi = np.arctan2(y, x)
            # This determines the profile we want for the skyrmion
            # Single twisting: k = pi / R
            k = pi_factor * np.pi / sk_r

            # We define here a 'hedgehog' skyrmion pointing down
            return (sign * np.sin(k * r) * np.cos(phi),
                    sign * np.sin(k * r) * np.sin(phi),
                    core * np.cos(k * r))
        else:
            if not out_skyrmion_dir:
                return (0, 0, -core)
            else:
                return out_skyrmion_dir

    # For TESTing purposes
    # This uses polygon mesh tools
    #def mu_s_in_hexagon(self, pos):
    #    x, y = pos[0], pos[1]

    #    if pmt.in_poly(x=x, y=y, n=6,
    #                   r=self.sim.mesh.Lx * 0.5,
    #                   translate=(self.sim.mesh.Lx * 0.5,
    #                              self.sim.mesh.Ly * 0.5)
    #                   ):

    #        return self.mu_s
    #    else:
    #        return 0

    # TODO: Compute the True helicoid wave length ?
    def helicoid_m_field(self, pos, _lambda=3.):
        """
        Generates a helicoid along the x direction with a default
        period of 3 nm (chirality according to the DMI)
        """
        x, y = pos[0], pos[1]

        return (-np.sin(np.pi * x / _lambda), 0, -np.cos(np.pi * x / _lambda))

    def sk_helicoid_m_field(self, pos, _lambda=2):
        """
        Generates a skyrmion and a helicoid in the same sample, using
        a period of 2 for the helicoid and the skyrmion located
        at 0.3 times the mesh size along x and 0.5 times along y
        """
        x, y = pos[0], pos[1]

        if x > 0.55 * self.sim.mesh.Lx:
            return self.helicoid_m_field(pos, _lambda=_lambda)
        else:
            return self.skyrmion_m_field(pos, sign=1,
                                         sk_pos=(0.3 * self.sim.mesh.Lx,
                                                 0.5 * self.sim.mesh.Ly)
                                         )

    def generate_simulation(self, do_precession=False, gamma=1.76e11,
                            load_mu_s=False):
        """

        Generates a simulation according to the self.simulation option

        gamma       :: Default is the free electron gyromagnetic ratio
        load_mu_s   :: A file path to a NPY file with the mu_s information
                       (only for the experimental_sample option)

        """

        # Mesh ----------------------------------------------------------------

        if self.simulation == 'experimental_sample':
            self.sim_from_image = sfi.sim_from_image(sfi.default_image)

            if not load_mu_s:
                self.sim_from_image.generate_magnetic_moments(mu_s=self.mu_s)
            else:
                self.sim_from_image.generate_magnetic_moments(
                    load_file=load_mu_s)

            self.sim = self.sim_from_image.sim

        elif self.simulation == '2D_square':
            # A square sized hexagonal mesh
            mesh = HexagonalMesh(self.mesh_a * 0.5, self.mesh_nx, self.mesh_ny,
                                 # periodicity=(True, True),
                                 alignment='square',
                                 unit_length=1e-9
                                 )
            self.sim = Sim(mesh)

            # If we use polygon mesh tools, we can use a hexagon shaped mesh
            # self.sim.mu_s = self.mu_s_in_hexagon

            self.sim.mu_s = self.mu_s

        elif self.simulation == '1D_chain':
            # A 1D chain using a cuboid mesh
            mesh = CuboidMesh(dx=self.mesh_a, nx=self.mesh_nx,
                              ny=1, nz=1,
                              # periodicity=(True, True),
                              unit_length=1e-9
                              )
            self.sim = Sim(mesh)

            # If we use polygon mesh tools, we can use a hexagon shaped mesh
            # self.sim.mu_s = self.mu_s_in_hexagon

            self.sim.mu_s = self.mu_s

        self.sim.driver.do_precession = do_precession
        self.sim.driver.gamma = gamma

        # Interactions --------------------------------------------------------

        exch = UniformExchange(self.J)
        self.sim.add(exch)

        dmi = DMI(D=(self.D), dmi_type='interfacial')
        self.sim.add(dmi)

        zeeman = Zeeman((self.B[0], self.B[1], self.B[2]))
        self.sim.add(zeeman, save_field=True)

        if self.ku:
            # Uniaxial anisotropy along + z-axis
            self.sim.add(Anisotropy(self.ku, axis=[0, 0, 1]))

        if self.Demag:
            print('Using Demag!')
            self.sim.add(DemagHexagonal())

        # ---------------------------------------------------------------------

        self.hls = np.ones_like(self.sim.spin.reshape(-1, 3))
        self.rgbs = np.ones((self.sim.spin.reshape(-1, 3).shape[0], 4))

    def update_DMI(self, D):
        self.D = D * const.meV
        self.sim.get_interaction('DMI').D = D * const.meV
        self.sim.driver.reset_integrator()

    def update_anisotropy(self, Ku):
        self.ku = Ku * const.meV
        self.sim.get_interaction('Anisotropy')._Ku = fidimag.common.helper.init_scalar(Ku * const.meV, self.sim.mesh)
        self.sim.driver.reset_integrator()

    def pin_boundaries(self, pin_direction=(0, 0, -1),
                       plot_component=None):
        """
        Pin the spins at the boundaries, setting the spin directions
        to *pin_direction*
        """
        boundary_spins = np.logical_and(np.any(self.sim.mesh.neighbours < 0,
                                               axis=1),
                                        self.sim.mu_s != 0)

        new_m = np.copy(self.sim.spin.reshape(-1, 3))
        new_m[boundary_spins] = np.array([pin_direction[0],
                                          pin_direction[1],
                                          pin_direction[2]
                                          ])
        self.sim.set_m(new_m.reshape(-1,))

        # Now we pin the spins:

        ngbs_filter = np.zeros(self.sim.pins.shape[0])
        # Filter rows by checking if any of the elements is less than zero
        # This means that if any of the neighbours of the i-th lattice site is
        # -1, we pin the spin direction at that site
        ngbs_filter = np.any(self.sim.mesh.neighbours < 0, axis=1,
                             out=ngbs_filter)

        self.sim.set_pins(ngbs_filter)

        # Plot the changes
        if plot_component:
            self.update_plot_component(plot_component)

    def release_boundaries(self):
        """
        Unpin the boundaries
        """
        ngbs_filter = self.sim.mu_s == 0
        self.sim.set_pins(ngbs_filter)

    def generate_m_field(self, m_function=None):
        """

        Requires that self.sim is defined
        This function generates the vector field for the magnetic moments

        If no function is provided, it automatically generates a skyrmion
        with the core pointing up, at the middle of the sample, for
        an interfacial system, where D > 0

        """
        if not m_function:
            m_function = lambda pos: self.skyrmion_m_field(pos, sign=1)

        self.sim.set_m(m_function)

    def plot_canvas(self, figsize=(8, 7)):
        """

        The base plot being updated during the simulation For the 2D
        simulations, we use a scatter plot. For the 1D chain we use a quiver
        plot showing the XZ plane. This plot is saved in the self.plot variable

        Spins are filtered to remove entities with zero magnetic moment mu_s

        """

        # An elongated figure for the 1D example
        if self.simulation == '1D_chain':
            self.fig = plt.figure(figsize=(8, 3))
        else:
            self.fig = plt.figure(figsize=figsize)

        self.ax = self.fig.add_subplot(111)

        # Assuming we have a homogeneous material, we scale the magnetic
        # moments before filtering since the magnitude is propertional to the
        # Bohr magneton, which is too small
        self._filter = self.sim.mu_s / self.mu_s > 1e-5

        if self.simulation == '1D_chain':
            self.plot = self.ax.quiver(
                # The chain is along X, Y is fixed
                self.sim.mesh.coordinates[:, 0][self._filter],  # X
                self.sim.mesh.coordinates[:, 1][self._filter],  # Y
                # Show the m_z components
                self.sim.spin.reshape(-1, 3)[:, 0][self._filter],  # m_x
                self.sim.spin.reshape(-1, 3)[:, 2][self._filter],  # m_z
                # color according to mz
                self.sim.spin.reshape(-1, 3)[:, 2][self._filter],
                cmap=mpl.cm.RdYlBu,
                pivot='mid', scale=1/0.04, linewidth=0.15
                )

            self.ax.set_xlabel(r'$x\,\, \mathrm{[nm]}$', fontsize=18)
            self.ax.set_ylabel(r'$z\,\, \mathrm{[nm]}$', fontsize=18)

            self.fig.subplots_adjust(top=0.85, bottom=0.2)

        else:
            self.plot = self.ax.scatter(
                self.sim.mesh.coordinates[:, 0][self._filter],
                self.sim.mesh.coordinates[:, 1][self._filter],
                # color according to mz
                c=self.sim.spin.reshape(-1, 3)[:, 2][self._filter],
                # Use hexagons as markers. The size will depend on the mesh
                # sizes
                s=40, marker='h', lw=0,
                cmap=mpl.cm.RdYlBu,
                vmin=-1, vmax=1,
                )
            self.ax.set_xlabel(r'$x\,\, \mathrm{[nm]}$', fontsize=18)
            self.ax.set_ylabel(r'$y\,\, \mathrm{[nm]}$', fontsize=18)

            self.fig.subplots_adjust(bottom=0.1)

        # Labels with simulation infor to update
        self.energy_text = self.ax.text(0, 1.05, '',
                                        transform=self.ax.transAxes,
                                        va='center', ha='left')

        self.step_text = self.ax.text(1., 1.05, '',
                                      transform=self.ax.transAxes,
                                      # Vertical and horizontal alignment
                                      va='center', ha='right')

        self.step_text_2 = self.ax.text(1., 1.08, '',
                                        transform=self.ax.transAxes,
                                        # Vertical and horizontal alignment
                                        va='bottom', ha='right')

        self.title_text = self.ax.text(0.5, 1.05, '',
                                       transform=self.ax.transAxes,
                                       # Vertical and horizontal alignment
                                       va='center', ha='center')

        # Set the ranges manually  if we use an external image to generate the
        # mesh
        if self.simulation == 'experimental_sample':
            # self.ax.imshow(self.sim_from_image.image_data,
            #                extent=self.sim_from_image.image_range)

            self.ax.set_xlim(self.sim_from_image.image_range[0],
                             self.sim_from_image.image_range[1]
                             )
            self.ax.set_ylim(self.sim_from_image.image_range[2],
                             self.sim_from_image.image_range[3]
                             )

        # Color bar -----------------------------------------------------------
        # Colour bar (append to not distort the main plot)
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)

        norm = mpl.colors.Normalize(vmin=-1, vmax=1)
        self.cbar = plt.colorbar(self.plot, cax=cax,
                                 cmap='RdYlBu',
                                 norm=norm,
                                 ticks=[-1, 0, 1],
                                 orientation='vertical',
                                 )
        self.cbar.set_label(r'$m_i$', rotation=270, labelpad=10, fontsize=16)

        # Interactive mode (not sure if necessary)
        plt.ion()

    def generate_widgets(self):
        """

        Generate the plot and the widgets with different options for
        the simulation

        """

        self.plot_canvas()

        # Select which m component to plot ------------------------------------

        mcomp_dropdown = ipw.Dropdown(options=['x', 'y', 'z', 'hls'],
                                      # description=r'$\text{Plot}\,\,m_{i}$',
                                      value='z'
                                      )
        mcomp_dropdown.observe(lambda b: self.static_plot_component_colorbar(
            m_component=mcomp_dropdown.value))
        mcomp_dropdown_text = ipw.Label(r'$\text{Plot}\,\,m_{i}$')
        mcomp_dropdown_text.width = '100px'

        mcomp_box = ipw.HBox(children=[mcomp_dropdown_text, mcomp_dropdown])

        # Button to start relaxation ------------------------------------------

        # We use two boxes with the maximum time and the number of steps
        run_box = ipw.FloatText()
        run_box_2 = ipw.FloatText()
        # Default values for the boxes
        run_box.value = 10
        run_box_2.value = 100
        # Change left and right margins
        run_box.layout.margin = '0px 50px'
        run_box.layout.width = '100px'
        run_box_2.layout.margin = '0px 50px'
        run_box_2.layout.width = '100px'

        # The button to start the relaxation of the system, using the values in
        # the boxes
        run_button = ipw.Button(description='Run sim!')
        run_button.layout.margin = '0px 50px'
        # The action requires a variable for the button (?). We just
        # impose it in a lambda function
        run_button.on_click(lambda b: self.relax(b,
                                                 max_time=run_box.value * 1e-9,
                                                 time_steps=run_box_2.value,
                                                 m_component=mcomp_dropdown.value
                                                 )
                            )
        # run_button.margin = '0px 20px'
        run_button.layout.background_color = '#9c0909'
        run_button.color = 'white'

        run_box_text = ipw.Label(r'$\text{Time [ns]}$')
        run_box_text.width = '50px'
        run_box_2_text = ipw.Label(r'$\text{Steps}$')
        run_box_2_text.width = '50px'

        # Align everything in a horizontal box
        run_container = ipw.HBox(children=[run_box_text, run_box,
                                           run_box_2_text, run_box_2,
                                           run_button]
                                 )

        # Boxes to update the Zeeman field ------------------------------------

        zeeman_box_texts = [ipw.Label(r'$B_x\,\,\text{[T]}$'),
                            ipw.Label(r'$B_y\,\,\text{[T]}$'),
                            ipw.Label(r'$B_z\,\,\text{[T]}$')
                            ]
        for text in zeeman_box_texts:
            text.width = '50px'

        zeeman_box_x = ipw.FloatText()
        zeeman_box_y = ipw.FloatText()
        zeeman_box_z = ipw.FloatText()

        # Show the default value in the box
        zeeman_box_x.value = self.B[0]
        zeeman_box_x.layout.width = '50px'
        zeeman_box_x.layout.margin = '0px 50px'
        zeeman_box_y.value = self.B[1]
        zeeman_box_y.layout.width = '50px'
        zeeman_box_y.layout.margin = '0px 50px'
        zeeman_box_z.value = self.B[2]
        zeeman_box_z.layout.width = '50px'
        zeeman_box_z.layout.margin = '0px 50px'

        # Update the simulation using a button
        zeeman_button = ipw.Button(description='Update field')
        zeeman_button.on_click(
            lambda button: self.update_Zeeman_field((zeeman_box_x.value,
                                                     zeeman_box_y.value,
                                                     zeeman_box_z.value,
                                                     )
                                                    )
            )

        # Draw the default Field in the plot title
        self.title_text.set_text('Field: {} T'.format(
            self.sim.get_interaction('Zeeman').B0))
        self.fig.canvas.draw()

        zeeman_container = ipw.HBox(children=[zeeman_box_texts[0],
                                              zeeman_box_x,
                                              zeeman_box_texts[1],
                                              zeeman_box_y,
                                              zeeman_box_texts[2],
                                              zeeman_box_z,
                                              zeeman_button
                                              ]
                                    )

        # DMI magnitude box ---------------------------------------------------

        DMI_box_text = ipw.Label(r'$\text{DMI [meV]}$')
        DMI_box_text.width = '50px'

        DMI_box = ipw.FloatText()

        # Show the default value in the box
        DMI_box.value = round(self.D / const.meV, 2)
        DMI_box.layout.width = '60px'
        DMI_box.layout.margin = '0px 50px'

        # Update the simulation using a button
        DMI_button = ipw.Button(description='Update')
        DMI_button.on_click(
            lambda button: self.update_DMI(DMI_box.value)
            )

        DMI_container = ipw.HBox(children=[DMI_box_text, DMI_box,
                                           DMI_button],
                                 )

        # Anisotropy -----------------------------------------------------------

        # anisotropy_box_text = ipw.Label(r'$\text{Anisotropy [meV]}$')
        # anisotropy_box_text.width = '100px'

        # anisotropy_box = ipw.FloatText()

        # # Show the default value in the box
        # anisotropy_box.value = round(self.ku / const.meV, 2)
        # anisotropy_box.layout.width = '60px'
        # anisotropy_box.layout.margin = '0px 50px'

        # # Update the simulation using a button
        # anisotropy_button = ipw.Button(description='Update')
        # anisotropy_button.on_click(
        #     lambda button: self.update_anisotropy(anisotropy_box.value)
        #     )

        # anisotropy_container = ipw.HBox(children=[anisotropy_box_text,
        #                                           anisotropy_box,
        #                                           anisotropy_button]
        #                                 )

        # Initial state -------------------------------------------------------

        # Options for the initial states. The values of the keys are the
        # initial magnetisation functions for Fidimag
        init_state_select = ipw.RadioButtons(
            options={'Skyrmion':
                     lambda pos: self.skyrmion_m_field(pos, sign=1),
                     '2-PI-Vortex':
                     lambda pos: self.skyrmion_m_field(pos, sign=1,
                                                       core=-1,
                                                       pi_factor=2,
                                                       out_skyrmion_dir=(0, 0, -1)
                                                       ),
                     '3-PI-Vortex':
                     lambda pos: self.skyrmion_m_field(pos, sign=1,
                                                       pi_factor=3,
                                                       ),
                     'Helicoid':
                     lambda pos: self.helicoid_m_field(pos),
                     'Sk-Helicoid':
                     lambda pos: self.sk_helicoid_m_field(pos),
                     'Random':
                     lambda pos: np.random.uniform(-1, 1, 3),
                     'Uniform':
                     (0, 0, -1)
                     },
            # description=r'$\text{Initial State}$'
            )
        init_state_select.selected_label = 'Skyrmion'

        # We need the extra variable for the buttons action
        def update_state(b):
            self.generate_m_field(m_function=init_state_select.value)
            self.update_plot_component(mcomp_dropdown.value)

        # The selection changes are taken with the observe method
        init_state_select.observe(update_state)

        init_state_select_text = ipw.Label(r'$\text{Initial State}$')
        init_state_select_text.width = '100px'

        init_state_select_box = ipw.HBox(children=[init_state_select_text,
                                                   init_state_select]
                                         )

        # Pin boundaries button -----------------------------------------------

        pin_button = ipw.Button(description='Pin boundaries')
        pin_button.on_click(
            lambda button: self.pin_boundaries(plot_component=mcomp_dropdown.value)
            )

        unpin_button = ipw.Button(description='Unpin boundaries')
        unpin_button.on_click(
            lambda button: self.release_boundaries()
            )

        pin_box = ipw.HBox(children=[pin_button, unpin_button])
        pin_box.layout.margin = '20px 0px'

        # Display -------------------------------------------------------------

        display(mcomp_box)
        display(zeeman_container)
        display(DMI_container)
        # display(anisotropy_container)
        display(init_state_select_box)
        display(pin_box)
        display(run_container)

        # return

    def relax(self, button,
              max_time=2e-9, time_steps=100,
              rtol=1e-8, atol=1e-10,
              m_component='z'
              ):
        """

        Relax the system until *max_time*, which has arbitrary units
        (we need to check this in Fidimag)

        time_steps indicate the number of steps to update the plot between 0
        and max_time. We will avoid updating the plot faster than 25 FPS, but
        we have seen that the delay is mostly from the run_until function, not
        from the plot redrawing process.  Thus decreasing the number of steps
        is more effective

        Tolerances are not being used for now

        """

        times = np.linspace(0, max_time, time_steps)

        # We won't call this function more than once every 1/25 seconds
        @throttle(seconds=1/25.)
        def local_plot_update(time):
            self.update_plot_component(m_component=m_component)
            # Update title
            self.step_text.set_text('Time: {:.5f} [ns]'.format(time * 1e9))

            # Show maximum dm/dt which seems more useful than the time
            # Avoid dividing by zero on the first steps
            if self.sim.driver.integrator.get_current_step() > 1e-12:
                self.step_text_2.set_text('Max dm/dt: {:.5f} [deg / ns]'.format(
                    self.sim.driver.compute_dmdt(self.sim.driver.integrator.get_current_step()) / self.DEGREE_PER_NANOSECOND
                    )
                    )
            # Update the energy
            self.energy_text.set_text('Energy: {:.4f} eV'.format(
                self.sim.compute_energy() / const.eV))

            # Update plot
            self.fig.canvas.draw()

        # Run the simulation --------------------------------------------------
        for time in times:

            self.sim.driver.run_until(time)

            # Update scatter or quiver plot
            local_plot_update(time)

        # Show the last update
        self.update_plot_component(m_component=m_component)

        # self.sim.save_vtk()

        self.sim.driver.reset_integrator()

    def update_Zeeman_field(self, B=(0, 0, 0)):
        self.B = B
        self.sim.get_interaction('Zeeman').update_field((B[0], B[1], B[2]))
        self.sim.driver.reset_integrator()

        # Update title on plot
        self.title_text.set_text('Field: {} T'.format(
            self.sim.get_interaction('Zeeman').B0)
            )
        self.fig.canvas.draw()

    def convert_to_rgb(self, hls_color):
        return np.array(colorsys.hls_to_rgb(hls_color[0] / (2 * np.pi),
                                            hls_color[1],
                                            hls_color[2]
                                            )
                        )

    def update_hls(self):
        self.hls[:, 0] = np.arctan2(self.sim.spin.reshape(-1, 3)[:, 1],
                                    self.sim.spin.reshape(-1, 3)[:, 0]
                                    )
        self.hls[:, 0][self.hls[:, 0] < 0] = self.hls[:, 0][self.hls[:, 0] < 0] + 2 * np.pi

        self.hls[:, 1] = 0.5 * (self.sim.spin.reshape(-1, 3)[:, 2] + 1)
        # self.hls[:, 2] = 0.5 * (self.sim.spin.reshape(-1, 3)[:, 2] + 1)

        self.rgbs = np.apply_along_axis(self.convert_to_rgb, 1, self.hls)
        # Some RGB values can get very small magnitudes, like 1e-10:
        self.rgbs[self.rgbs < 0] = 0

    def update_plot_component(self, m_component='z'):
        """
        Update the vector data for the plot (the spins do not move
        so we don't need to update the coordinates) and redraw
        """
        m = self.sim.spin.reshape(-1, 3)
        if self.simulation == '1D_chain':
            self.plot.set_UVC(m[:, 0], m[:, 2], m[:, 2])
        else:
            if m_component == 'hls':
                self.update_hls()
                # self.plot.set_array(self.hls[:, 0][self._filter])
                self.plot.set_color(self.rgbs[self._filter])
            else:
                self.plot.set_array(m[:, self.v_dict[m_component]][self._filter])
        self.fig.canvas.draw()

    def static_plot_component_colorbar(self, m_component='z'):
        if m_component == 'hls':
            self.plot.set_cmap('hsv')
            self.cbar.set_ticklabels([r'$0$', '', r'$2\pi$'])
            self.fig.canvas.draw()
        else:
            self.cbar.set_ticklabels(['-1', '0', '1'])
            self.plot.set_cmap('RdYlBu')
            self.fig.canvas.draw()

        self.update_plot_component(m_component)

    def relax_sim(self):
        self.sim.driver.do_precession = False
        self.sim.driver.relax(dt=1e-13)
