/*
 * Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
 *
 * Copyright 2020- The Mem3DG Authors
 * and the project initiators Cuncheng Zhu, Christopher T. Lee, and
 * Padmini Rangamani.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Please help us support Mem3DG development by citing the research
 * papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
 * for more information.
 */

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
                                                      Settings for mesh mutator
                                                     )delim");
  meshmutator.def(py::init<>(),
                  R"delim(
                    Constructor of mesh mutator object
                  )delim");
  meshmutator.def_readwrite("mutateMeshPeriod",
                            &MeshProcessor::MeshMutator::mutateMeshPeriod,
                            R"delim(
        Frequency (in numbers of iterations) at which to apply of mesh mutations.
      )delim");
  meshmutator.def_readonly("isFlipEdge",
                           &MeshProcessor::MeshMutator::isFlipEdge,
                           R"delim(
                             Whether to perform edge flipping
                           )delim");
  meshmutator.def_readonly("isSplitEdge",
                           &MeshProcessor::MeshMutator::isSplitEdge,
                           R"delim(
                              Whether edges should be split.

                              Note:
                                  Increases the number of elements in the mesh
                           )delim");
  meshmutator.def_readonly("isCollapseEdge",
                           &MeshProcessor::MeshMutator::isCollapseEdge,
                           R"delim(
        Whether edges should be collapsed.

        Note:
            Decreases the number of elements in the mesh
      )delim");
  meshmutator.def_readonly("isChangeTopology",
                           &MeshProcessor::MeshMutator::isChangeTopology,
                           R"delim(
                              Whether the mesh manipulations changing topology are permitted
                           )delim");

  /**
   * @brief vertex shifting
   */
  meshmutator.def_readwrite("isShiftVertex",
                            &MeshProcessor::MeshMutator::isShiftVertex,
                            R"delim(
                              Boolean controlling if vertex shifting is run
                            )delim");

  /**
   * @brief mesh smoothing
   */
  meshmutator.def_readwrite("isSmoothenMesh",
                            &MeshProcessor::MeshMutator::isSmoothenMesh,
                            R"delim(
                              Boolean controlling if mesh smoothing operations are run
                            )delim");

  /**
   * @brief flipping criterion
   */
  meshmutator.def_readwrite("flipNonDelaunay",
                            &MeshProcessor::MeshMutator::flipNonDelaunay,
                            R"delim(
                              Whether to flip non - Delaunay edge
                            )delim");
  meshmutator.def_readwrite(
      "flipNonDelaunayRequireFlat",
      &MeshProcessor::MeshMutator::flipNonDelaunayRequireFlat,
      R"delim(
        Whether to require flatness condition when flipping non - Delaunay edge
      )delim");

  /**
   * @brief splitting criterion
   */
  meshmutator.def_readwrite("splitLarge",
                            &MeshProcessor::MeshMutator::splitLarge,
                            R"delim(
                              Split edges with large incident faces
                            )delim");
  meshmutator.def_readwrite("splitLong", &MeshProcessor::MeshMutator::splitLong,
                            R"delim(
                              Split long edge
                            )delim");
  meshmutator.def_readwrite("splitCurved",
                            &MeshProcessor::MeshMutator::splitCurved,
                            R"delim(
                              Split edges on high curvature domains
                            )delim");
  meshmutator.def_readwrite("splitCurvedScaleFactor",
                            &MeshProcessor::MeshMutator::splitCurvedScaleFactor,
                            R"delim(
        Split edge if it is longer than this scale factor times the estimated optimal edge length given by curvTol.
    )delim");

  meshmutator.def_readwrite("splitSharp",
                            &MeshProcessor::MeshMutator::splitSharp,
                            R"delim(
                              Split edge with sharp membrane property changes

                              .. warning::
                                  Current not implemented
                            )delim");
  meshmutator.def_readwrite("splitFat", &MeshProcessor::MeshMutator::splitFat,
                            R"delim(
                              Split obtuse triangles
                            )delim");
  meshmutator.def_readwrite("splitSkinnyDelaunay",
                            &MeshProcessor::MeshMutator::splitSkinnyDelaunay,
                            R"delim(
        Split triangles with poor aspect ratio which are still Delaunay
      )delim");
  meshmutator.def_readwrite("minimumEdgeLength",
                            &MeshProcessor::MeshMutator::minimumEdgeLength,
                            R"delim(
                              Target minimum edge length
                            )delim");
  meshmutator.def_readwrite("maximumEdgeLength",
                            &MeshProcessor::MeshMutator::maximumEdgeLength,
                            R"delim(
                              Target maximum edge length
                            )delim");

  /**
   * @brief collapsing criterion
   */
  meshmutator.def_readwrite("collapseSkinny",
                            &MeshProcessor::MeshMutator::collapseSkinny,
                            R"delim(
                              Collapse skinny triangles
                            )delim");
  meshmutator.def_readwrite("collapseSmall",
                            &MeshProcessor::MeshMutator::collapseSmall,
                            R"delim(
                              Collapse small triangles
                            )delim");
  meshmutator.def_readwrite("collapseFlat",
                            &MeshProcessor::MeshMutator::collapseFlat,
                            R"delim(
                              Collapse flat edges
                            )delim");
  meshmutator.def_readwrite(
      "collapseFlatScaleFactor",
      &MeshProcessor::MeshMutator::collapseFlatScaleFactor,
      R"delim(
        Collapse edge if it is shorter than this scale factor times the
            estimated optimal edge length given by curvTol.
      )delim");
  meshmutator.def_readwrite("curvTol", &MeshProcessor::MeshMutator::curvTol,
                            R"delim(
                              Tolerance for curvature approximation
                            )delim");

  meshmutator.def_readwrite("maximumFaceArea",
                            &MeshProcessor::MeshMutator::maximumFaceArea,
                            R"delim(
                              Target maximum face area
                            )delim");

  meshmutator.def_readwrite("minimumFaceArea",
                            &MeshProcessor::MeshMutator::minimumFaceArea,
                            R"delim(
                              Target minimum face area
                            )delim");

  py::class_<MeshProcessor> meshprocessor(pymem3dg, "MeshProcessor",
                                          R"delim(
                                            Object to store mesh processing settings
                                          )delim");
  meshprocessor.def(py::init<>(),
                    R"delim(
                      MeshProcessor constructor
                    )delim");
  meshprocessor.def_readwrite("meshMutator", &MeshProcessor::meshMutator,
                              R"delim(
                                meshMutator object
                              )delim");
  meshprocessor.def_readonly("isMeshMutate", &MeshProcessor::isMeshMutate,
                             R"delim(
                               Whether to perform mesh mutation operations
                             )delim");
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
