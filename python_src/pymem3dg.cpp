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
#include <cstdarg>
#include <cstddef>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Eigen/src/Core/util/Constants.h"

#include "mem3dg/solver/geometry.h"
#include "mem3dg/solver/integrator/conjugate_gradient.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/version.h"

#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include "mem3dg/mem3dg"
#include "pybind11/cast.h"

#ifdef MEM3DG_WITH_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

namespace py = pybind11;
namespace gc = ::geometrycentral;

namespace mem3dg {
namespace solver {
namespace integrator {

// Forward declarations
void init_integrators(py::module_ &);
void init_forces(py::module_ &);
void init_mutators(py::module_ &);
void init_geometry(py::module_ &);
void init_system(py::module_ &);
void init_parameters(py::module_ &);
void init_energy(py::module_ &);

// Initialize the `pymem3dg` module
PYBIND11_MODULE(_core, pymem3dg) {
  pymem3dg.doc() = "Python wrapper around the DDG solver C++ library.";
  pymem3dg.attr("__version__") = MEM3DG_VERSION;
  pymem3dg.attr("__short_version__") = MEM3DG_VERSION_SHORT;

  init_integrators(pymem3dg);
  init_forces(pymem3dg);
  init_mutators(pymem3dg);
  init_geometry(pymem3dg);
  init_system(pymem3dg);
  init_parameters(pymem3dg);
  init_energy(pymem3dg);

#ifdef MEM3DG_WITH_GPERFTOOLS
  pymem3dg.def("startProfiler", &ProfilerStart, "start profiler",
               py::arg("filename"));
  pymem3dg.def("stopProfiler", &ProfilerStop, "stop profiler");
#else
  pymem3dg.def(
      "startProfiler",
      []([[maybe_unused]] std::string foo) {
        mem3dg_runtime_warning("Library is not linked with gperftools. "
                               "Profiling is not available.");
      },
      "start profiler", py::arg("filename"));
  pymem3dg.def(
      "stopProfiler",
      []() {
        mem3dg_runtime_warning("Library is not linked with gperftools. "
                               "Profiling is not available.");
      },
      "stop profiler");
#endif

#pragma region mesh_io
  // ==========================================================
  // =============      mesh generation    ===============
  // ==========================================================
  pymem3dg.def("getCylinder", &getCylinderMatrix,
               "get face and vertex matrix of a non-capped cylinder",
               py::arg("radius"), py::arg("radialSubdivision"),
               py::arg("axialSubdivision"), py::arg("frequency") = 1,
               py::arg("amplitude") = 0);

  pymem3dg.def("getIcosphere", &getIcosphereMatrix,
               "get face and vertex matrix of icosphere", py::arg("radius"),
               py::arg("subdivision") = 0);

  pymem3dg.def("getTetrahedron", &getTetrahedronMatrix,
               "get face and vertex matrix of tetrahedron");

  pymem3dg.def("getDiamond", &getDiamondMatrix,
               "get face and vertex matrix of diamond", py::arg("dihedral"));

  pymem3dg.def("getHexagon", &getHexagonMatrix,
               "get face and vertex matrix of Hexagon", py::arg("radius"),
               py::arg("subdivision") = 0);

  pymem3dg.def(
      "linearSubdivide",
      py::overload_cast<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &,
                        Eigen::Matrix<double, Eigen::Dynamic, 3> &,
                        std::size_t>(&linearSubdivide),
      py::arg("face"), py::arg("vertex"), py::arg("nSub"), R"delim(
          subdivide the mesh with linear interpolation
      )delim");

  pymem3dg.def(
      "loopSubdivide",
      py::overload_cast<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &,
                        Eigen::Matrix<double, Eigen::Dynamic, 3> &,
                        std::size_t>(&loopSubdivide),
      "subdivide the mesh in Loop scheme", py::arg("face"), py::arg("vertex"),
      py::arg("nSub"), R"delim(
          subdivide the mesh in loop scheme
      )delim");

  pymem3dg.def("getFaceAndVertexMatrix", &getFaceAndVertexMatrix,
               py::arg("plyName"), R"delim(
          read face matrix and vertex matrix from .ply file
      )delim");
  pymem3dg.def(
      "getDataElementName", &getDataElementName, py::arg("plyName"),
      R"delim(retrieve all richData element name from .ply file. Namely the list of the places where data live in, such as vertex, edge or face.
      )delim");
  pymem3dg.def(
      "getDataPropertyName", &getDataPropertyName, py::arg("plyName"),
      py::arg("elementName"),
      R"delim(retrieve all richData property name from .ply file. Namely the list of the data where data on the particular element, such as vertex, edge or face.
      )delim");
  pymem3dg.def("getData", &getData, py::arg("plyName"), py::arg("elementName"),
               py::arg("propertyName"),
               R"delim(read richData from .ply file
      )delim");
  pymem3dg.def(
      "processSoup", &processSoup, py::arg("meshName"),
      R"delim( process soup data in .ply and return face and vertex matrices.
      )delim");
  pymem3dg.def(
      "getFaceSurfacePointClosestToEmbeddedCoordinate",
      py::overload_cast<const EigenVectorX3sr &, const EigenVectorX3dr &,
                        const std::array<double, 3> &,
                        const Eigen::Matrix<bool, Eigen::Dynamic, 1> &,
                        const std::array<bool, 3> &>(
          &getFaceSurfacePointClosestToEmbeddedCoordinate),
      py::arg("faceMatrix"), py::arg("vertexMatrix"),
      py::arg("embeddedCoordinate"), py::arg("filter"),
      py::arg("accountedCoordinate") = std::array<bool, 3>{true, true, true},
      R"delim(find the face surface point closest to a embedded coordinate)delim");
  //   pymem3dg.def(
  //       "getFaceSurfacePointClosestToEmbeddedCoordinate",
  //       py::overload_cast<const EigenVectorX3sr &, const EigenVectorX3dr &,
  //                         const std::array<double, 3> &,
  //                         const std::array<bool, 3> &>(
  //           &getFaceSurfacePointClosestToEmbeddedCoordinate),
  //       py::arg("faceMatrix"), py::arg("vertexMatrix"),
  //       py::arg("embeddedCoordinate"),
  //       py::arg("accountedCoordinate") = std::array<bool, 3>{true, true,
  //       true}, R"delim(find the face surface point closest to a embedded
  //       coordinate)delim");
  pymem3dg.def(
      "getVertexClosestToEmbeddedCoordinate",
      py::overload_cast<const EigenVectorX3dr &, const std::array<double, 3> &,
                        const Eigen::Matrix<bool, Eigen::Dynamic, 1> &,
                        const std::array<bool, 3> &>(
          &getVertexClosestToEmbeddedCoordinate),
      py::arg("vertexMatrix"), py::arg("embeddedCoordinate"), py::arg("filter"),
      py::arg("accountedCoordinate") = std::array<bool, 3>{true, true, true},
      R"delim(find the index of vertex closest to a embedded coordinate)delim");
  //   pymem3dg.def(
  //       "getVertexClosestToEmbeddedCoordinate",
  //       py::overload_cast<const EigenVectorX3dr &, const std::array<double,
  //       3> &,
  //                         const std::array<bool, 3> &>(
  //           &getVertexClosestToEmbeddedCoordinate),
  //       py::arg("vertexMatrix"), py::arg("embeddedCoordinate"),
  //       py::arg("accountedCoordinate") = std::array<bool, 3>{true, true,
  //       true}, R"delim(find the index of vertex closest to a embedded
  //       coordinate)delim");
  pymem3dg.def(
      "getVertexFurthestFromBoundary",
      py::overload_cast<const EigenVectorX3sr &, const EigenVectorX3dr &>(
          &getVertexFurthestFromBoundary),
      py::arg("faceMatrix"), py::arg("vertexMatrix"),
      R"delim(find the vertex furthest away from the boundaries)delim");

#pragma endregion mesh_io
};
} // namespace integrator
} // namespace solver
} // namespace mem3dg
