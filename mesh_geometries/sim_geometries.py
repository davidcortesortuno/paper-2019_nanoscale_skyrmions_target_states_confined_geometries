from __future__ import print_function

import fidimag
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
import numpy as np
from matplotlib.path import Path
import fidimag.common.constant as const
# import polygonmeshtools as pmt
import matplotlib.pyplot as plt
import matplotlib


class TruncatedTriangleSim(object):
    """
    Class to generate a simulation object on a mesh with a truncated
    triangle geometry.

    To generate the material inside the triangular geometry, we use
    the Path class from matplotlib, which has a .contains_point method
    to check if points are inside irregular polygons
    """

    def __init__(self, L, offsets, lattice_a=0.2715, mu_s=3,
                 name='truncated_triangle_sim', shells=1, driver='llg'):
        self.L = L
        self.shells = shells

        if isinstance(offsets, float) or isinstance(offsets, int):
            self.offsets = [offsets, offsets, offsets]
        else:
            self.offsets = offsets

        self.lattice_a = lattice_a
        self.mu_s = mu_s * const.mu_B
        self.name = name
        self.driver = driver

        # Generate the simulation with a truncated triangle as a mesh
        self.generate_triangle()
        self.generate_mesh()
        self.generate_simulation()

    def generate_triangle(self):
        """
        Generates a Matplotlib Path object which is a polygon defined by
        a series of points
                                                      _
                              --             .         \
                             |              ...         \  offset_1
                             |         C   .....   D    _\
                             |            o*****o
                             |           *********
          L * sin(pi / 3)    |          ***********
                             |         *************
                             |    B  o***************o  E
                             |      ..***************..
                             |     ....o***********o....
                             --    |   A           F
                                 (0,0)

                                  |____|           |____|
                                  offset_0         offset_2

                                  |_____________________|
                                             L
        """

        codes = []
        vertices = []

        # Vertices:

        # Shift the x coordinates by half the amount of the offset to the left,
        # so we avoid creating a lot of points without material at the left
        # side of the triangle
        xshift = 0.5 * self.offsets[0] - self.lattice_a
        # xshift = 0

        # A
        vertices.append([self.offsets[0] - xshift,
                         0])
        # B
        vertices.append([self.offsets[0] * 0.5 - xshift,
                         self.offsets[0] * np.sin(np.pi / 3)])
        # C
        vertices.append([self.L * 0.5 - self.offsets[1] * 0.5 - xshift,
                         self.L * np.sin(np.pi / 3)
                         - self.offsets[1] * np.sin(np.pi / 3)
                         ])
        # D
        vertices.append([self.L * 0.5 + self.offsets[1] * 0.5 - xshift,
                         self.L * np.sin(np.pi / 3)
                         - self.offsets[1] * np.sin(np.pi / 3)
                         ])
        # E
        vertices.append([self.L - self.offsets[2] * 0.5 - xshift,
                         self.offsets[2] * np.sin(np.pi / 3)
                         ])
        # F
        vertices.append([self.L - self.offsets[2] - xshift,
                         0])

        vertices.append([self.offsets[0] - xshift,
                         0])  # Ignored

        # Codes:
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]

        self.vertices, self.codes = vertices, codes
        self.truncated_triangle = Path(self.vertices, self.codes)

    def generate_mesh(self):

        # The width of the mesh will be the width of the triangle less
        # the offsets (that's why we shifted the material to the left)
        # We add 2 to avoid cropping the mesh if the integer approximation
        # is not optimal
        nx = int(np.around((self.L - (self.offsets[0] + self.offsets[2]) * 0.5)
                           / self.lattice_a)) + 2

        # Height of the triangle divided by the lattice hexagons height
        # +1 is to avoid cropping the triangle
        ny = int(np.around((self.L - self.offsets[1]) * np.sin(np.pi / 3)
                           / (self.lattice_a * np.sqrt(3) * 0.5))
                 ) + 1

        self.mesh = HexagonalMesh(self.lattice_a * 0.5, nx, ny,
                                  alignment='square', unit_length=1e-9,
                                  shells=self.shells
                                  )

    def inside_triangle(self, r):
        """
        This function requires that the Path self.truncated_triangle
        is already defined
        """

        if self.truncated_triangle.contains_point([r[0], r[1]]):
            return self.mu_s
        else:
            return 0

    def generate_simulation(self):
        """
        This requires that the mesh is defined
        """
        self.sim = fidimag.atomistic.Sim(self.mesh, name=self.name,
                                         driver=self.driver)
        self.sim.set_mu_s(self.inside_triangle)
        self.sim.set_m((0, 0, 1))

    def plot_triangle_matplotlib(self):
        patch = matplotlib.patches.PathPatch(self.truncated_triangle,
                                             facecolor='orange', lw=2)
        f = plt.figure(figsize=(12, 12))
        ax = f.add_subplot(111)
        ax.add_patch(patch)
        ax.set_xlim(0, self.L)
        ax.set_ylim(0, self.L)
        plt.show()


class HexagonSim(object):
    """
    Class to generate a simulation object on a hexagonally shaped mesh.

    To generate the material inside the hexagon, we use the polygon-mesh-tools
    from Github
    """

    def __init__(self, R, lattice_a=0.2715, mu_s=3,
                 name='hexagon_sim', shells=1, driver='llg'):

        self.shells = shells
        self.R = R
        self.lattice_a = lattice_a
        self.mu_s = mu_s * const.mu_B
        self.name = name
        self.driver = driver

        # Generate the simulation with a truncated triangle as a mesh
        self.generate_hexagon()
        self.generate_mesh()
        self.generate_simulation()

    def generate_hexagon(self):
        """
        Generates a Matplotlib Path object which is a polygon defined by
        a series of points
        """

        codes = []
        vertices = []

        # Vertices:

        # A
        vertices.append([0,
                         self.R * np.cos(np.pi / 6)])
        # B
        vertices.append([self.R * (1 - np.sin(np.pi / 6)),
                         2 * self.R * np.cos(np.pi / 6)])
        # C
        vertices.append([self.R * (1 + np.sin(np.pi / 6)),
                         2 * self.R * np.cos(np.pi / 6)])
        # D
        vertices.append([2 * self.R,
                         self.R * np.cos(np.pi / 6)])
        # E
        vertices.append([self.R * (1 + np.sin(np.pi / 6)),
                         0])
        # F
        vertices.append([self.R * (1 - np.sin(np.pi / 6)),
                         0])

        vertices.append([0,
                         self.R * np.cos(np.pi / 6)])  # Ignored

        # Codes:
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]

        self.vertices, self.codes = vertices, codes
        self.hexagon = Path(self.vertices, self.codes)

    def generate_mesh(self):

        nx = int(np.around(2 * self.R / self.lattice_a)) + 1
        ny = int(np.around(2 * self.R * np.cos(np.pi / 6)
                           / (self.lattice_a * np.sqrt(3) * 0.5)
                           )) + 1

        self.mesh = HexagonalMesh(self.lattice_a * 0.5, nx, ny,
                                  alignment='square', unit_length=1e-9,
                                  shells=self.shells
                                  )

    # def inside_hexagon(self, r):
    #     """
    #     """
    #     if pmt.in_poly(x=r[0], y=r[1], n=6,
    #                    r=self.R,
    #                    translate=(self.mesh.Lx * 0.5,
    #                               self.mesh.Ly * 0.5
    #                               )
    #                    ):

    #         return self.mu_s
    #     else:
    #         return 0

    def inside_hexagon(self, r):
        """
        """
        if self.hexagon.contains_point([r[0], r[1]]):
            return self.mu_s
        else:
            return 0

    def generate_simulation(self):
        """
        This requires that the mesh is defined
        """
        self.sim = fidimag.atomistic.Sim(self.mesh, name=self.name,
                                         driver=self.driver)
        self.sim.set_mu_s(self.inside_hexagon)
        self.sim.set_m((0, 0, 1))

if __name__ == '__main__':
    # An example (test?)
    triang = TruncatedTriangleSim(15, 2)
    triang.sim.save_vtk()

    triang = TruncatedTriangleSim(15, [3, 3, 5],
                                  name='truncated_triangle_assymetric')
    triang.sim.save_vtk()

    hexagon = HexagonSim(10)
    hexagon.sim.save_vtk()
