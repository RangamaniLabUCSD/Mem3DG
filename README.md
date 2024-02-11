# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)

[![Testing and release](https://github.com/RangamaniLabUCSD/Mem3DG/actions/workflows/ci.yaml/badge.svg?branch=main)](https://github.com/RangamaniLabUCSD/Mem3DG/actions/workflows/ci.yaml)
[![Zenodo](https://zenodo.org/badge/244037679.svg)](https://zenodo.org/doi/10.5281/zenodo.10359392)
[![PyPI](https://img.shields.io/pypi/v/pymem3dg)](https://pypi.org/project/pymem3dg/)
[![Conda-forge](https://anaconda.org/conda-forge/pymem3dg/badges/version.svg)](https://anaconda.org/conda-forge/pymem3dg)

Mem3DG is a flexible software package to model the membrane and its dynamics using unstructured meshes.
This work is currently under heavy development, please star this repository to follow along!

## Installation

The core of Mem3DG is written in C++ with functions exposed to a Python interface using pybind11.
The python interface and helper utilities for visualization and analysis are bundled together into a python package called pymem3dg.
pymem3dg can be obtained using pip (`pip install pymem3dg`) or conda (`conda install pymem3dg`) from their respective repositories.
For the majority of users who wish to develop models, this is the recommended installation procedure.

### Dependencies

We acknowledge the use of helpful external libraries including:

* [Geometry-Central](https://geometry-central.net/), the core mesh data structure
* [libigl](https://libigl.github.io/), some geometry input generation
* [pybind11](https://pybind11.readthedocs.io/en/stable/), for Python-C++ interoperability
* [Polyscope](https://polyscope.run/py/), trajectory visualization GUI
* [PCG](https://www.pcg-random.org/index.html), random number generation
* [NetCDF-cxx4](https://github.com/Unidata/netcdf-cxx4), binary trajectory I/O

While the majority of these libraries are either included in this repository

### Building for development

For advanced users who seek to modify the code, develop simulations in C++ directly, or produce/link to Mem3DG shared libraries, Mem3DG can be configured and built using CMake.
After cloning the repository, initialize submodules and their dependencies, and follow the standard CMake out-of-source configuration and build procedures.

```bash
git submodule update --init --recursive
mkdir build
cd build
cmake -DWITH_NETCDF=ON -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release
```

There are several CMake options which you can configure.
Further options can be discovered by inspecting the CMakeLists.txt or using a utility like `ccmake` or `cmake -LAH`.
These options can be passed during the configuration step (e.g., `-DWITH_NETCDF=ON`).
```cmake
option(BUILD_PYMEM3DG "Build the python extensions?" ON)
option(WITH_NETCDF "Build with NetCDF (binary trajectory output)?" ON)
option(BUILD_MEM3DG_DOCS "Configure documentation" OFF)
option(M3DG_GET_OWN_EIGEN "Download own Eigen" ON)
option(M3DG_GET_OWN_PYBIND11 "Download own pybind11" ON)
option(WITH_LIBUNWIND "Link libunwind for stack traces" OFF)
option(WITH_OPENMP "Build with OpenMP support" OFF)
option(LINK_PROFILER "Link profiler gperftools" OFF)
option(MEM3DG_PEDANTIC "Be extremely pedantic while compiling" OFF)
```

### Building `pymem3dg`

While the pymem3dg extension module can be built directly with CMake, packaging and installing in a matter which conforms with guidelines from the [Python packaging authority](https://packaging.python.org/en/latest/) is best done using python tooling.
We use [`scikit-build-core`](https://scikit-build-core.readthedocs.io/en/latest/) to bridge between CMake and the Python build system.
Typical Python metadata and details are specified in `pyproject.toml`.

pymem3dg can be built and installed by using pip from the root of this repository.
```bash
git submodule update --init --recursive
pip install -v .
```
Extra dependencies for building the documentation and tests can be installed as well `pip install -v .[docs,tests]`.
Other options can be passed to CMake by modifying config-settings in pip.
`pip install . --config-settings=cmake.build-type="Debug"`

## Acknowledging your use of Mem3DG

Mem3DG is developed by [Cuncheng Zhu](https://github.com/cuzhucuncheng), [Christopher T. Lee](https://ctlee.github.io/), with contributions from others.
Development of Mem3DG is funded in part by AFOSR MURI FA9550-18-1-0051, and a Hartwell Foundation Postdoctoral Fellowship.
