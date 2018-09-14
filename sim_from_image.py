from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh

import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

default_image = ('Romming_data/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw/Pd_Fe_Ir111_20160118_CrTip_015_B-4T_003_rotated_large_bw.png')


class sim_from_image(object):
    def __init__(self,
                 image_path,
                 image_range=[0, 21.03, 0, 17.79],
                 mesh_a=0.2715,
                 shells=1,
                 # mesh_nx=79,
                 # mesh_ny=78
                 sim_name='unnamed',
                 bg_colour=(1., 1., 1.),
                 driver='llg'
                 ):

        """

        Generate an object on top of the fidimag simulation class, to generate
        an atomistic simulation using a hexagonal mesh, where the magnetic
        moments are populated in the mesh according to the pixels in a PNG
        image, located in *image_path* The mesh size is generated automatically
        from the dimensions provided for the image (the *image_range* array),
        so that the mesh covers the whole picture. It is recommended that the
        image has white borders since mesh points away from the picture range
        will take the values from the closest boundary point (see the
        populate_mesh function)

        image_range  :: An array [xmin, xmax, ymin, ymax] with the ranges
                        of the image, given in nm

        """

        # print 'Mesh spacings are: nx = {}, ny = {}'.format(mesh_a,
        #                                                    mesh_a * np.sqrt(3) * 0.5)

        self.image_path = image_path
        self.image_range = image_range

        mesh_nx = int(np.around(np.abs(image_range[1] - image_range[0]) / mesh_a))
        mesh_ny = int(np.around(
            np.abs(image_range[3] - image_range[2]) / (mesh_a * np.sqrt(3) * 0.5)
            ))

        self.mesh = HexagonalMesh(mesh_a * 0.5, mesh_nx, mesh_ny,
                                  # periodicity=(True, True),
                                  alignment='square',
                                  unit_length=1e-9,
                                  shells=shells
                                  )

        self.sim = Sim(self.mesh, name=sim_name, driver=driver)

        # Generate the 2D matrix with the image
        # If the image is B/W, every entry will have either [1., 1., 1.]
        # or [0., 0., 0.] representing the RGB data from every pixel
        self.image_data = mpl.image.imread(self.image_path)

        # We will create two new matrices, with the X and Y coordinates
        # of every point in the image matrix, using the range
        # provided in the arguments. For doing this, we use a meshgrid, which
        # automatically generates a 2D mesh of coordinates. The range
        # in the Y coordinates matrix is reversed since the image matrix
        # starts from top to bottom to generate the picture.
        # The shape is [n_rows, n_columns]  --> (ny, nx)
        self.image_coords_x = np.linspace(self.image_range[0],
                                          self.image_range[1],
                                          self.image_data.shape[1])

        self.image_coords_y = np.linspace(self.image_range[3],
                                          self.image_range[2],
                                          self.image_data.shape[0])

        (self.image_coords_x,
         self.image_coords_y) = np.meshgrid(self.image_coords_x,
                                            self.image_coords_y)

        self.bg_colour = bg_colour

    def populate_mesh(self, pos, mu_s):
        """
        Magnetisation scalar field function to fill the hexagonal mesh
        according to the position of the pixels in the image matrix.

        For a given point of the mesh, this function computes the
        minimum distance to a point in the image matrix, whose coordinates
        are given by the mesh matrices image_coords_x and image_coords_y
        Then it finds the position in the matrix (i, j),
        where the closest point lies and check if that pixel is black
        (for now we assume any other pixel is white, i.e. [1., 1., 1.]).
        If so, it returns a magnetic moment value that point

        Very far away points from the range of the image, the points in
        the mesh will take the values from the boundary points. Thus
        it is recommended to use an image with a white border (?)

        """
        x, y = pos[0], pos[1]

        # Distance in X and Y (Euclidean)
        dist_x = self.image_coords_x - x
        dist_y = self.image_coords_y - y
        # A matrix with the distances of this point (x, y) with respect
        # to every point in the image matrix coordinates
        dist = np.sqrt(dist_x ** 2 + dist_y ** 2)

        i, j = np.where(dist == np.min(dist))
        i = i[0]
        j = j[0]

        # Set material for any colour that is not the background colour
        if not (abs(self.image_data[i, j][0] - self.bg_colour[0]) < 1e-5 and
                abs(self.image_data[i, j][1] - self.bg_colour[1]) < 1e-5 and
                abs(self.image_data[i, j][2] - self.bg_colour[2]) < 1e-5):
            return mu_s
        else:
            return 0

    def generate_magnetic_moments(self, mu_s=1,
                                  load_file=False,
                                  save_file=False
                                  ):
        """
        save    :: A file name to save the magnetic moments array
        load    :: A path to a file with the array of magnetic moments
                   (obtained with the *save* option)
        """

        if not load_file:
            self.sim.set_mu_s(lambda pos: self.populate_mesh(pos, mu_s))
        else:
            self.sim.set_mu_s(np.load(load_file))

        if save_file:
            np.save(save_file, self.sim.mu_s)

    def plot_mesh(self, figsize=(20, 7.5)):
        """
        Plot the mesh on top of the image
        """
        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_subplot(111)

        ax1.scatter(self.mesh.coordinates[:, 0],
                    self.mesh.coordinates[:, 1],
                    c=self.sim.mu_s,
                    s=35,
                    marker='h', lw=0,
                    cmap=mpl.cm.Reds,
                    # vmin=0, vmax=self.sim.mu_s,
                    # gridsize=(63, 41)
                    )

        ax1.imshow(self.image_data, extent=self.image_range)

        ax1.set_xlim([self.image_range[0],
                      self.image_range[1]]
                     )
        ax1.set_ylim([self.image_range[2],
                      self.image_range[3]]
                     )

    # .........................................................................

    def populate_spin_orientations(self, pos, orientations, colours):
        """
        """
        x, y = pos[0], pos[1]

        # Distance in X and Y (Euclidean)
        dist_x = self.image_coords_x - x
        dist_y = self.image_coords_y - y
        # A matrix with the distances of this point (x, y) with respect
        # to every point in the image matrix coordinates
        dist = np.sqrt(dist_x ** 2 + dist_y ** 2)

        i, j = np.where(dist == np.min(dist))
        i = i[0]
        j = j[0]

        for k, spin_direction in enumerate(orientations):
            if (abs(self.image_data[i, j][0] - colours[k][0]) < 1e-5 and
                abs(self.image_data[i, j][1] - colours[k][1]) < 1e-5 and
                abs(self.image_data[i, j][2] - colours[k][2]) < 1e-5):

                return spin_direction

        return (0., 0., 0.)

    def generate_magnetic_configuration(self, orientations, colours):
        """
        """

        self.sim.set_m(lambda pos: self.populate_spin_orientations(pos,
                                                                   orientations,
                                                                   colours))
