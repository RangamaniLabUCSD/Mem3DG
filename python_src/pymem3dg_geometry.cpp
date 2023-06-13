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

#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include "mem3dg/mem3dg"
#include "pybind11/cast.h"

namespace gc = ::geometrycentral;
namespace mem3dg {
namespace solver {
namespace integrator {
namespace py = pybind11;

void init_geometry(py::module_ &pymem3dg) {
  // ==========================================================
  // =============          Geometry            ===============
  // ==========================================================
  py::class_<Geometry> geometry(pymem3dg, "Geometry",
                                R"delim(
        Geometry class
    )delim");

  /**
   * @brief Constructors
   */
  geometry.def(
      py::init<std::string, Eigen::Matrix<bool, Eigen::Dynamic, 1> &>(),
      py::arg("inputMesh"), py::arg("notableVertex"),
      R"delim(
        Geometry constructor with .ply files.
      )delim");
  geometry.def(py::init<std::string>(), py::arg("inputMesh"),
               R"delim(
          Geometry constructor with .ply files.
        )delim");
  geometry.def(py::init<EigenVectorX3sr &, EigenVectorX3dr &,
                        Eigen::Matrix<bool, Eigen::Dynamic, 1> &>(),
               py::arg("faceMatrix"), py::arg("vertexMatrix"),
               py::arg("notableVertex"),

               R"delim(
        Geometry constructor with Matrices
      )delim");
  geometry.def(py::init<EigenVectorX3sr &, EigenVectorX3dr &>(),
               py::arg("faceMatrix"), py::arg("vertexMatrix"),
               R"delim(
        Geometry constructor with Matrices
      )delim");
#ifdef MEM3DG_WITH_NETCDF
  geometry.def(py::init<std::string, int>(), py::arg("trajFile"),
               py::arg("startingFrame"),
               R"delim(
        Geometry constructor with NetCDF trajectory file
      )delim");
#endif
  /**
   * @brief Getters
   */
  geometry.def(
      "getSurfaceArea", [](Geometry &s) { return s.surfaceArea; },
      py::return_value_policy::copy,
      R"delim(
          get the surface area of the mesh
      )delim");
  geometry.def(
      "getGeodesicDistance",
      [](Geometry &s) { return s.geodesicDistance.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the geodesic distance from notable vertices of the mesh
      )delim");
  geometry.def(
      "getVolume", [](Geometry &s) { return s.volume; },
      py::return_value_policy::copy,
      R"delim(
          get the enclosed volume of the mesh
      )delim");
  geometry.def(
      "getNotableVertex", [](Geometry &s) { return s.notableVertex.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get vertex data to track the notable vertex, which may or may not be a single vertex
      )delim");
  geometry.def(
      "getLumpedMassMatrix",
      [](Geometry &s) { return s.vpg->vertexLumpedMassMatrix; },
      py::return_value_policy::copy,
      R"delim(
          get the lumped mass matrix of the mesh
      )delim");
  geometry.def(
      "getCotanLaplacian", [](Geometry &s) { return s.vpg->cotanLaplacian; },
      py::return_value_policy::copy,
      R"delim(
          get the Cotan Laplacian matrix of the mesh
      )delim");
  geometry.def(
      "getVertexNormals",
      [](Geometry &s) { return toMatrix(s.vpg->vertexNormals); },
      py::return_value_policy::copy,
      R"delim(
          get angle-weighted normal on vertices
      )delim");
  geometry.def(
      "getVertexMatrix",
      [](Geometry &s) {
        return gc::EigenMap<double, 3>(s.vpg->inputVertexPositions);
      },
      py::return_value_policy::copy,
      R"delim(
          get the vertex matrix
      )delim");
  geometry.def(
      "getFaceMatrix",
      [](Geometry &s) { return s.mesh->getFaceVertexMatrix<std::size_t>(); },
      py::return_value_policy::copy,
      R"delim(
          get the face matrix
      )delim");
  geometry.def(
      "getVertexAdjacencyMatrix", [](Geometry &s) { return s.vpg->d0; },
      py::return_value_policy::copy,
      R"delim(
          get the signed E-V vertex adjacency matrix, equivalent of d0 operator
      )delim");
  geometry.def(
      "getEdgeAdjacencyMatrix", [](Geometry &s) { return s.vpg->d1; },
      py::return_value_policy::copy,
      R"delim(
          get the signed F-E edge adjacency matrix, equivalent of d1 operator
      )delim");
  geometry.def(
      "getEdgeLengths", [](Geometry &s) { return s.vpg->edgeLengths.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get edge lengths
      )delim");
  geometry.def(
      "getVertexDualAreas",
      [](Geometry &s) { return s.vpg->vertexDualAreas.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get vertex dual area
      )delim");
  geometry.def(
      "getFaceAreas", [](Geometry &s) { return s.vpg->faceAreas.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get face area
      )delim");
  geometry.def(
      "getVertexMeanCurvatures",
      [](Geometry &s) {
        s.vpg->requireVertexMeanCurvatures();
        return s.vpg->vertexMeanCurvatures.raw();
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated scalar mean curvature
      )delim");
  geometry.def(
      "getVertexGaussianCurvatures",
      [](Geometry &s) {
        s.vpg->requireVertexGaussianCurvatures();
        return s.vpg->vertexGaussianCurvatures.raw();
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated scalar Gaussian Curvature
      )delim");
  geometry.def(
      "getVertexGaussianCurvatureVectors",
      [](Geometry &s) {
        auto vector = s.computeVertexGaussianCurvatureVectors();
        return toMatrix(vector);
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated vector Gaussian Curvature
      )delim");
  geometry.def(
      "getVertexVolumeVariationVectors",
      [](Geometry &s) {
        auto vector = s.computeVertexVolumeVariationVectors();
        return toMatrix(vector);
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated vector volume variation (dual area)
      )delim");
  geometry.def(
      "getVertexMeanCurvatureVectors",
      [](Geometry &s) {
        auto vector = s.computeVertexMeanCurvatureVectors();
        return toMatrix(vector);
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated vector Mean Curvature
      )delim");
  geometry.def(
      "getVertexSchlafliLaplacianMeanCurvatureVectors",
      [](Geometry &s, EigenVectorX1d spontaneousCurvature) {
        gcs::VertexData<double> H0(*s.mesh);
        H0.fromVector(spontaneousCurvature);
        auto vector = s.computeVertexSchlafliLaplacianMeanCurvatureVectors(H0);
        return toMatrix(vector);
      },
      py::return_value_policy::copy,
      R"delim(
          get the vertex Schlafli based Laplacian of mean curvature Vectors
      )delim");
  geometry.def(
      "getPolyscopePermutations",
      [](Geometry &s) { return gcs::polyscopePermutations(*s.mesh); },
      R"delim(
          get polyscope permutation
      )delim");
  geometry.def(
      "getPolyscopeEdgeOrientations",
      [](Geometry &s) { return gcs::polyscopeEdgeOrientations(*s.mesh).raw(); },
      R"delim(
          get polyscope edge orientation
      )delim");

  /**
   * @brief    geometry setter
   */
  geometry.def(
      "setInputVertexPositions",
      [](Geometry &s, EigenVectorX3dr newGeo) {
        gc::EigenMap<double, 3>(s.vpg->inputVertexPositions) = newGeo;
      },
      R"delim(
          set the vertex matrix
      )delim");

  /**
   * @brief Methods
   */
  geometry.def("computeGeodesicDistance", &Geometry::computeGeodesicDistance,
               R"delim(

          compute the geodesic distance centered around Center cached in System

          return: NDarray[double]
      )delim");
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
