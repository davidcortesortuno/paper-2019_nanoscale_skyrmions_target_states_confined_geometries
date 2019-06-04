FROM ubuntu:18.04

# -----------------------------------------------------------------------------
# LIBRARIES

# ENTRYPOINT ["/bin/bash"]  
# SHELL ["/bin/bash", "-c"] 

# Avoid user interaction dialog
ENV DEBIAN_FRONTEND=noninteractive

# We need cython3 from the ubuntu repos, since cython-ising Sundials fail
# with pip cython
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt install -y build-essential cmake python3-dev python3-pip cython3 \
    python3-pytest-cov liblapack-dev libopenblas-dev \
    wget make gfortran libblas-dev liblapack-dev python3-tk sudo fonts-lato

RUN ln -s /usr/bin/python3 /usr/bin/python

# Next line for codecov target
RUN apt install -y curl git

RUN pip3 install ipywidgets nbval numpy scipy matplotlib notebook psutil pytest pytest-cov pyvtk -U
# RUN pip3 install --upgrade setuptools==20.4

# Headless Matplotlib:
ENV MPLBACKEND=Agg

# FIDIMAG ---------------------------------------------------------------------

WORKDIR /usr/local
RUN git clone https://github.com/computationalmodelling/fidimag.git
WORKDIR /usr/local/fidimag
# Work with stable release
RUN git checkout tags/v2.9

# Install CVODE and FFTW libraries
WORKDIR /usr/local/fidimag/bin
RUN bash install-fftw.sh
RUN bash install-sundials.sh

ENV PYTHONPATH="/usr/local/fidimag:$PYTHONPATH"
ENV LD_LIBRARY_PATH="/usr/local/fidimag/local/lib:$LD_LIBRARY_PATH"

WORKDIR /usr/local/fidimag
RUN python3 setup.py build_ext --inplace

# -----------------------------------------------------------------------------

# Set threads for OpenMP:
ENV OMP_NUM_THREADS=2
# WORKDIR /io

# -----------------------------------------------------------------------------
# User to make Binder happy
ENV NB_USER micromag
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

# Set the working directory
WORKDIR /home/${USER}/micromag

# -----------------------------------------------------------------------------
