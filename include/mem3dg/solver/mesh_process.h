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

struct MeshProcessor {
  struct MeshRegularizer {
    /// Whether regularize mesh
    bool isMeshRegularize = false;

    /// triangle ratio constant
    double Kst = 0;
    /// Local stretching modulus
    double Ksl = 0;
    /// Edge spring constant
    double Kse = 0;
    /// whether vertex shift
    bool shiftVertex = true;

    /**
     * @brief summarizeStatus
     */
    void summarizeStatus() {
      isMeshRegularize = (Kst != 0) || (Ksl != 0) || (Kse != 0);
    };
  };

  struct MeshMutator {
    /// Whether mutate mesh
    bool isMeshMutate = false;

    /// Whether edge flip
    bool isEdgeFlip = false;
    /// Whether split edge
    bool isSplitEdge = false;
    /// Whether collapse edge
    bool isCollapseEdge = false;

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
    /// split poor aspected triangle that is still Delaunay
    bool splitSkinnyDelaunay = false;

    /// collapse skinny triangles
    bool collapseSkinny = false;
    /// collapse small triangles
    bool collapseSmall = false;
    /// whether require flatness condition when collapsing small edge
    bool collapseSmallNeedFlat = false;

    /// tolerance for curvature approximation
    double curvTol = 0.0012;
    /// target face area
    double targetFaceArea = 0.001;

    /**
     * @brief summarizeStatus
     */
    void summarizeStatus() {
      isEdgeFlip = (flipNonDelaunay || flipNonDelaunayRequireFlat);
      isSplitEdge = (splitCurved || splitLarge || splitLong || splitSharp ||
                     splitSkinnyDelaunay);
      isCollapseEdge =
          (collapseSkinny || collapseSmall || collapseSmallNeedFlat);
      isMeshMutate = isEdgeFlip || isSplitEdge || isCollapseEdge;
    };

    /**
     * @brief return condition of edge flip
     */
    bool ifFlip(const gcs::Edge e, const gcs::VertexPositionGeometry &vpg);

    /**
     * @brief return condition of edge split
     */
    bool ifSplit(const gc::Edge e, const gcs::VertexPositionGeometry &vpg);

    /**
     * @brief return condition of edge collapse
     */
    bool ifCollapse(const gc::Edge e, const gcs::VertexPositionGeometry &vpg);

    void maskAllNeighboring(gcs::VertexData<bool> &smoothingMask,
                            const gcs::Vertex v);

    void neighborAreaSum(const gcs::Edge e,
                         const gcs::VertexPositionGeometry &vpg, double &area,
                         std::size_t &num_neighbor);

    double
    computeCurvatureThresholdLength(const gcs::Edge e,
                                    const gcs::VertexPositionGeometry &vpg);
  };

  MeshMutator meshMutator;
  MeshRegularizer meshRegularizer;

  /**
   * @brief summarizeStatus
   */
  void summarizeStatus() {
    meshMutator.summarizeStatus();
    meshRegularizer.summarizeStatus();
    std::cout << " and " << (meshMutator.isMeshMutate && meshRegularizer.isMeshRegularize) << std::endl;
    if (meshMutator.isMeshMutate && meshRegularizer.isMeshRegularize) {
      mem3dg_runtime_error("For topology changing simulation, mesh "
                           "regularization cannot be applied!");
    }
  };
};

} // namespace solver
} // namespace mem3dg
