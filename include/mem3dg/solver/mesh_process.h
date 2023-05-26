// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2021:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#pragma once

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include <math.h>
#include <vector>

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"
#include "mem3dg/constants.h"
#include "mem3dg/macros.h"
#include "mem3dg/mesh_io.h"
#include "mem3dg/meshops.h"
#include "mem3dg/type_utilities.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {
namespace solver {
/**
 * @brief Mesh processor
 *
 */
struct MeshProcessor {

  /**
   * @brief Mesh mutator object
   *
   */
  struct MeshMutator {
    /// period of mesh mutation
    std::size_t mutateMeshPeriod = std::numeric_limits<std::size_t>::max();

    /// Whether edge flip
    bool isFlipEdge = false;
    /// Whether split edge
    bool isSplitEdge = false;
    /// Whether collapse edge
    bool isCollapseEdge = false;
    /// Whether change topology
    bool isChangeTopology = false;

    /// whether vertex shift
    bool isShiftVertex = false;

    // whether conduct smoothing
    bool isSmoothenMesh = false;

    /// flip non-delaunay edge
    bool flipNonDelaunay = false;
    /// whether require flatness condition when flipping non-Delaunay edge
    bool flipNonDelaunayRequireFlat = false;

    /// split edge with large faces
    bool splitLarge = false;
    /// split long edge
    bool splitLong = false;
    /// split edge on high curvature domain
    bool splitCurved = false;
    /// split edge with sharp membrane property change
    bool splitSharp = false;
    /// split obtuse triangle
    bool splitFat = false;
    /// split poor aspect ratio triangle that is still Delaunay
    bool splitSkinnyDelaunay = false;
    /// min edge length
    double minimumEdgeLength = 0.001;
    /// max edge length
    double maximumEdgeLength = 0.01;

    /// collapse skinny triangles based on angles
    bool collapseSkinny = false;
    /// collapse triangles smaller than the \a targetFaceArea
    bool collapseSmall = false;
    /// target face area
    double targetFaceArea = 0.001;
    /// whether require flatness condition when collapsing small edge
    bool collapseFlat = false;

    /// tolerance for curvature approximation
    double curvTol = 0.0012;

    /**
     * @brief Aggregate fine mutation flags to set broader flag states
     */
    void summarizeStatus();

    /**
     * @brief Check flip conditions for edge e
     *
     * Affected by \a flipNonDelaunay and \a flipNonDelaunayRequireFlat flags
     *
     * @param e edge of interest
     * @param vpg geometry
     * @return true if edge \a e should be flipped
     * @return false if edge \a e is ok
     */
    bool checkFlipCondition(const gcs::Edge e,
                            const gcs::VertexPositionGeometry &vpg);

    /**
     * @brief return condition of edge split
     */
    bool checkSplitCondition(const gc::Edge e,
                             const gcs::VertexPositionGeometry &vpg);

    /**
     * @brief return condition of edge collapse
     */
    bool checkCollapseCondition(const gc::Edge e,
                                const gcs::VertexPositionGeometry &vpg);

    void markVertices(gcs::VertexData<bool> &mutationMarker,
                      const gcs::Vertex v, const size_t layer = 0);

    std::tuple<double, std::size_t>
    neighborAreaSum(const gcs::Edge e, const gcs::VertexPositionGeometry &vpg);

    double
    computeCurvatureThresholdLength(const gcs::Edge e,
                                    const gcs::VertexPositionGeometry &vpg);
  };

  /// mesh mutator
  MeshMutator meshMutator;
  /// Whether mutate mesh
  bool isMeshMutate = false;

  /**
   * @brief Aggregate fine mutation flags to set broader flag states
   */
  void summarizeStatus();
};

} // namespace solver
} // namespace mem3dg
