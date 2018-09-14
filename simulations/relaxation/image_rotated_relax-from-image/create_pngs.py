import os,shutil
import glob
import re

# First, and before importing any Enthought packages, set the ETS_TOOLKIT
# environment variable to qt4, to tell Traits that we will use Qt.
# os.environ['ETS_TOOLKIT'] = 'qt4'
# print os.getenv("ETS_TOOLKIT")

from mayavi import mlab
import mayavi
import matplotlib
from matplotlib.cm import datad, get_cmap
import numpy as np

vtk_folders_path = 'vtks/'

# List of vtk files
vtk_folders = glob.glob(vtk_folders_path + '*')

# imrootd = os.path.join('png', '{}'.format(args.sim_name))
imrootd = 'pngs/'

if not os.path.exists(imrootd):
    os.makedirs(imrootd)

# Figure specs
f = mlab.figure(bgcolor=(1, 1, 1),
                fgcolor=(0, 0, 0),
                size=(600, 600)
                )

f.scene.off_screen_rendering = True

for vtkf in vtk_folders:

    sim_name = vtkf[:-5]
    vtk_file = glob.glob(vtkf + '/*.vtk')[0]
    data = mlab.pipeline.open(vtk_file)

    # We will clasify the files by the radius R and initial state
    state = re.search('(?<=n-)[a-z-\d]+', sim_name).group(0)

    png_folder = 'pngs'
    if not os.path.exists(png_folder):
        os.makedirs(png_folder)

    # Extract vector norm and filter the points whose norm is
    # zero (we can set an arbitrary low value)
    vnorm = mlab.pipeline.extract_vector_norm(data)

    vtres = mlab.pipeline.threshold(vnorm)
    try:
        vtres.lower_threshold = 1e-5
    except:
        print('No points with zero vector norm')

    # Extract vec comp and plot
    vecomp = mlab.pipeline.extract_vector_components(vtres)

    # Extract z-component of the data
    vecomp.component = 'z-component'

    surf = mlab.pipeline.surface(vecomp,
                                 vmax=1, vmin=-1,
                                 colormap='gist_earth'
                                 )
    surf.actor.property.interpolation = 'flat'

    # if args.reversed_map:
    #     surf.module_manager.scalar_lut_manager.reverse_lut = True

    # View from top as default
    mlab.view(elevation=0,
              azimuth=0,
              # distance=0
              )
    # f.scene.z_plus_view()

    try:
        # If there is no point with zero 'm', this threshold fails
        vtres.lower_threshold = 1e-5
    except:
        pass

    # Zoom for saving purposes
    f.scene.magnification = 2

    # Save in an appropriate file

    f.scene.save(os.path.join(png_folder, '{}.png'.format(state)))

    mlab.clf(figure=f)
