// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// 
// Copyright (c) 2020:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
// 

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mem3dg/solver/ddgsolver.h"

#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>

#include "mem3dg/solver/force.h"

namespace ddgsolver {
namespace py = pybind11;

// Initialize the `pymem3dg` module
PYBIND11_MODULE(pymem3dg, pymem3dg) {
  pymem3dg.doc() = "Python wrapper around the DDG solver C++ library.";

  pymem3dg.def("driver_ply", &driver_ply, "the driver function for input .ply mesh file ",
            py::arg("inputMesh"), py::arg("refMesh"), py::arg("isTuftedLaplacian"), py::arg("isProtein"),
            py::arg("mollifyFactor"), py::arg("isVertexShift"), py::arg("Kb"), py::arg("H0"), py::arg("sharpness"),
            py::arg("r_H0"), py::arg("Kse"), py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"),
            py::arg("Kv"), py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"), py::arg("gamma"), 
            py::arg("kt"), py::arg("ptInd"), py::arg("Kf"), 
            py::arg("conc"), py::arg("height"), py::arg("radius"), py::arg("h"), py::arg("T"), 
            py::arg("eps"), py::arg("closeZone"), py::arg("increment"), 
            py::arg("tSave"), py::arg("tMollify"), py::arg("outputDir"),
            R"delim(
               Run the driver.

               Args:
                   filename (:py:class:`str`): Filename to read.
            
               Returns:
                   :py:class:`int`: success.
            )delim");

    pymem3dg.def("driver_nc", &driver_nc, " a driver function for input .nc trajectory file",
      py::arg("trajFile"), py::arg("startingFrame"), py::arg("isTuftedLaplacian"), py::arg("isProtein"),
      py::arg("mollifyFactor"), py::arg("isVertexShift"), py::arg("Kb"),
      py::arg("H0"), py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"),
      py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"),
      py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"), py::arg("gamma"),
      py::arg("kt"), py::arg("ptInd"), py::arg("Kf"), py::arg("conc"),
      py::arg("height"), py::arg("radius"), py::arg("h"), py::arg("T"),
      py::arg("eps"), py::arg("closeZone"), py::arg("increment"),
      py::arg("tSave"), py::arg("tMollify"), py::arg("outputDir"),
      R"delim(
               Run the driver.

               Args:
                   filename (:py:class:`str`): Filename to read.
            
               Returns:
                   :py:class:`int`: success.
            )delim");

  pymem3dg.def("viewer", &viewer, " a visualization function",
            py::arg("fileName"));


  pymem3dg.def("view_animation", &view_animation, " a visualization function",
            py::arg("fileName"));


  pymem3dg.def("genIcosphere", &genIcosphere, "Generate a icosphere .ply file",
            py::arg("nSub"), py::arg("path"), py::arg("R"));
};
} // namespace ddgsolver
