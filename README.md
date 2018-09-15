# Simulation of Pd/Fe/Ir(111) systems

![](images/hex_island_sk_helix.png)

This repository contains the libraries and scripts to completely reproduce the simulations of the publication *...*. These simulations are based on the [Fidimag](http://doi.org/10.5334/jors.223) code [2] which can run discrete spin simulations.

Images of the experimental data used in the simulations are located in the `Romming_data` folder.

## Scripts

Simulation scripts are located in the `simulation` folder, where they are separated in three categories:

- **relaxation**: scripts that directly relax a given initial state for a particular geometry (the experimental island, hexagons or truncated triangles). Magnetic parameters can also be modified, see the `bash` scripts for details. The library containing all the available options and initial states is `hexagonal_fidimag.py`, which uses `argparse`.  The *relaxation* directory also contains the scripts to fully reproduce the phase diagrams shown in [1]. These simulations require substantial simulation time and disk space, thus it is recommended to start them with precaution.

- **hysteresis**: these scripts take the initial state from the *relaxation* simulations (specific paths to `npy` files might require updating the final step number) and perform a sequential field sweep of the system, according to the specified options. Scripts are written as `cfg` Python config files and the main library with the options is `hexagonal_hysteresis.py`.

- **NEBM**: these scripts also rely on the simulations from the *relaxation* folder to specify the initial state. The main library in this case is `hexagonal_neb_fidimag.py`, which uses the GNEBM implementation from Fidimag.

A library to create a mesh/simulation from an image is given in `sim_from_image.py`. In our case we use the islands from experiments in the `Romming_data` folder. A copy of the numpy array containing the spin data, which was btained with this library for the island system, is stored in the `mu_s` directory.

A library to create meshes with different geometries is given in the `mesh_geometries` folder.



# Widget

An interactive widget is provided in the `simulation_widget_nb.ipynb` notebook which runs the `simulation_widget.py` library. This widget is a proof of concept for the extensibility of Fidimag, in this case with the IPython widgets library. The notebook shows a simulation for the island system reproduced from experiments but further options are available in `simulation_widget.py`

# References

[1]

[2] Bisotti, M.-A., Cortés-Ortuño, D., Pepper, R., Wang, W., Beg, M., Kluyver, T., & Fangohr, H. (2018). Fidimag – A Finite Difference Atomistic and Micromagnetic Simulation Package. Journal of Open Research Software, 6(1), 22. DOI: http://doi.org/10.5334/jors.223
