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

namespace gc = ::geometrycentral;
namespace mem3dg {
namespace solver {
namespace integrator {
namespace py = pybind11;

void init_mutators(py::module_ &pymem3dg) {
  // ==========================================================
  // =============     MeshMutator              ===============
  // ==========================================================
  py::class_<MeshProcessor::MeshMutator> meshmutator(pymem3dg, "MeshMutator",
                                                     R"delim(
        The mesh mutator settings
    )delim");
  meshmutator.def(py::init<>(),
                  R"delim(
       meshmutator constructor
      )delim");
  meshmutator.def_readwrite("mutateMeshPeriod",
                            &MeshProcessor::MeshMutator::mutateMeshPeriod,
                            R"delim(
          period of mesh mutation. measured in # of integration iteration
      )delim");
  meshmutator.def_readonly("isFlipEdge",
                           &MeshProcessor::MeshMutator::isFlipEdge,
                           R"delim(
          get the option of whether do edge flip
      )delim");
  meshmutator.def_readonly("isSplitEdge",
                           &MeshProcessor::MeshMutator::isSplitEdge,
                           R"delim(
          get the option of whether split edge to grow mesh
      )delim");
  meshmutator.def_readonly("isCollapseEdge",
                           &MeshProcessor::MeshMutator::isCollapseEdge,
                           R"delim(
          get the option of whether Collapse edge to grow mesh
      )delim");
  meshmutator.def_readonly("isChangeTopology",
                           &MeshProcessor::MeshMutator::isChangeTopology,
                           R"delim(
          get the option of change topology
      )delim");

  /**
   * @brief vertex shifting
   */
  meshmutator.def_readwrite("isShiftVertex",
                            &MeshProcessor::MeshMutator::isShiftVertex,
                            R"delim(
          get the option of whether do vertex shift
      )delim");

  /**
   * @brief mesh smoothing
   */
  meshmutator.def_readwrite("isSmoothenMesh",
                            &MeshProcessor::MeshMutator::isSmoothenMesh,
                            R"delim(
          whether conduct mesh smoothing operation
      )delim");

  /**
   * @brief flipping criterion
   */
  meshmutator.def_readwrite("flipNonDelaunay",
                            &MeshProcessor::MeshMutator::flipNonDelaunay,
                            R"delim(
          whether flip non-Delaunay edge
      )delim");
  meshmutator.def_readwrite(
      "flipNonDelaunayRequireFlat",
      &MeshProcessor::MeshMutator::flipNonDelaunayRequireFlat,
      R"delim(
          whether require flatness condition when flipping non-Delaunay edge
      )delim");

  /**
   * @brief splitting criterion
   */
  meshmutator.def_readwrite("splitLarge",
                            &MeshProcessor::MeshMutator::splitLarge,
                            R"delim(
          split edge with large faces
      )delim");
  meshmutator.def_readwrite("splitLong", &MeshProcessor::MeshMutator::splitLong,
                            R"delim(
          split long edge
      )delim");
  meshmutator.def_readwrite("splitCurved",
                            &MeshProcessor::MeshMutator::splitCurved,
                            R"delim(
          split edge on high curvature domain
      )delim");
  meshmutator.def_readwrite("splitSharp",
                            &MeshProcessor::MeshMutator::splitSharp,
                            R"delim(
          split edge with sharp membrane property change
      )delim");
  meshmutator.def_readwrite("splitFat", &MeshProcessor::MeshMutator::splitFat,
                            R"delim(
          split obtuse triangle
      )delim");
  meshmutator.def_readwrite("splitSkinnyDelaunay",
                            &MeshProcessor::MeshMutator::splitSkinnyDelaunay,
                            R"delim(
          split poor aspect ratio triangle that is still Delaunay
      )delim");
  meshmutator.def_readwrite("minimumEdgeLength",
                            &MeshProcessor::MeshMutator::minimumEdgeLength,
                            R"delim(
          minimum edge length
      )delim");
  meshmutator.def_readwrite("maximumEdgeLength",
                            &MeshProcessor::MeshMutator::maximumEdgeLength,
                            R"delim(
          maximum edge length
      )delim");

  /**
   * @brief collapsing criterion
   */
  meshmutator.def_readwrite("collapseSkinny",
                            &MeshProcessor::MeshMutator::collapseSkinny,
                            R"delim(
          collapse skinny triangles
      )delim");
  meshmutator.def_readwrite("collapseSmall",
                            &MeshProcessor::MeshMutator::collapseSmall,
                            R"delim(
          collapse small triangles
      )delim");
  meshmutator.def_readwrite("collapseFlat",
                            &MeshProcessor::MeshMutator::collapseFlat,
                            R"delim(
         collapse flat edge
      )delim");
  meshmutator.def_readwrite("curvTol", &MeshProcessor::MeshMutator::curvTol,
                            R"delim(
          tolerance for curvature approximation
      )delim");
  meshmutator.def_readwrite("targetFaceArea",
                            &MeshProcessor::MeshMutator::targetFaceArea,
                            R"delim(
          target face area
      )delim");

  py::class_<MeshProcessor> meshprocessor(pymem3dg, "MeshProcessor",
                                          R"delim(
        The mesh processor settings
    )delim");
  meshprocessor.def(py::init<>(),
                    R"delim(
       meshprocessor constructor
      )delim");
  meshprocessor.def_readwrite("meshMutator", &MeshProcessor::meshMutator,
                              R"delim(
          meshMutator struct
      )delim");
  meshprocessor.def_readonly("isMeshMutate", &MeshProcessor::isMeshMutate,
                             R"delim(
          get the option of whether do mesh mutation
      )delim");
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
