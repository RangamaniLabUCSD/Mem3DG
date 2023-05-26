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
#include "mem3dg/solver/mesh_process.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"
#include "geometrycentral/utilities/vector3.h"
#include "mem3dg/constants.h"
#include "mem3dg/meshops.h"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void MeshProcessor::summarizeStatus() {
  meshMutator.summarizeStatus();
  isMeshMutate = meshMutator.isChangeTopology || meshMutator.isShiftVertex ||
                 meshMutator.isSmoothenMesh;
};

bool MeshProcessor::MeshMutator::checkFlipCondition(
    const gcs::Edge e, const gcs::VertexPositionGeometry &vpg) {
  gcs::Halfedge he = e.halfedge();

  if (flipNonDelaunay && !e.isBoundary()) {
    // Check if flatness condition is required
    if (flipNonDelaunayRequireFlat) {
      if (!(abs(vpg.edgeDihedralAngle(he.edge())) < (constants::PI / 36))) {
        // short circuit if flatness required and not flat
        return false;
      }
    }
    // True if non-delaunay
    return (vpg.cornerAngle(he.next().next().corner()) +
            vpg.cornerAngle(he.twin().next().next().corner())) >
           (constants::PI);
  }
  return false;
}

void MeshProcessor::MeshMutator::summarizeStatus() {
  isFlipEdge = (flipNonDelaunay || flipNonDelaunayRequireFlat);
  isSplitEdge = (splitCurved || splitLarge || splitLong || splitSharp ||
                 splitSkinnyDelaunay);
  isCollapseEdge = (collapseSkinny || collapseSmall || collapseFlat);
  isChangeTopology = isFlipEdge || isSplitEdge || isCollapseEdge;
};

bool MeshProcessor::MeshMutator::checkCollapseCondition(
    const gc::Edge e, const gcs::VertexPositionGeometry &vpg) {
  if (vpg.edgeLength(e) >= maximumEdgeLength) {
    return false;
  }

  gcs::Halfedge he = e.halfedge();
  bool isBoundary = e.isBoundary();
  if (!he.isInterior()) {
    he = he.twin();
  }
  bool condition = false;

  bool is2Small = false;
  bool is2Skinny = false;
  bool isFlat = false;
  bool isSmooth = true;
  bool isCollapse = false;

  // conditions for collapsing
  if (collapseSkinny) {
    // is2Skinny = (vpg.cornerAngle(he.next().next().corner()) +
    //              vpg.cornerAngle(he.twin().next().next().corner())) <
    //             constants::PI / 3;

    if (isBoundary) {
      // check opposite angle
      is2Skinny =
          (vpg.cornerAngle(he.next().next().corner()) < constants::PI / 6);
    } else {
      is2Skinny =
          (vpg.cornerAngle(he.next().next().corner()) +
           vpg.cornerAngle(he.twin().next().next().corner())) <
              constants::PI / 3 ||
          (vpg.cornerAngle(he.next().corner()) +
           vpg.cornerAngle(he.twin().corner())) > (constants::PI * 1.333) ||
          (vpg.cornerAngle(he.corner()) +
           vpg.cornerAngle(he.twin().next().corner())) >
              (constants::PI * 1.333);
    }
    // if (is2Skinny) {
    //   mem3dg_print("is2Skinny");
    // }
    condition = is2Skinny;
  }

  if (collapseSmall) {
    double areaSum;
    std::size_t num_neighbor;
    std::tie(areaSum, num_neighbor) = neighborAreaSum(e, vpg);

    if (isBoundary) {
      is2Small = (areaSum - vpg.faceArea(he.face())) <
                 (num_neighbor - 1) * targetFaceArea;
    } else {
      is2Small =
          (areaSum - vpg.faceArea(he.face()) - vpg.faceArea(he.twin().face())) <
          (num_neighbor - 2) * targetFaceArea;
    }
    // isSmooth =
    //     abs(H0[he.tipVertex()] - H0[he.tailVertex()]) < (0.667 *
    //     targetdH0);
    if (is2Small) {
      if (abs(vpg.edgeDihedralAngle(e)) > (constants::PI / 6)) {
        is2Small = false;
      }
    }
    condition = condition || is2Small;
  }

  if (collapseFlat) {
    isFlat =
        vpg.edgeLength(e) < (0.667 * computeCurvatureThresholdLength(e, vpg));
    condition = condition || isFlat;
  }

  return condition;
}

bool MeshProcessor::MeshMutator::checkSplitCondition(
    const gcs::Edge e, const gcs::VertexPositionGeometry &vpg) {

  if (vpg.edgeLength(e) <= minimumEdgeLength) {
    return false;
  }

  gcs::Halfedge he = e.halfedge();
  bool isBoundary = e.isBoundary();
  if (!he.isInterior()) {
    he = he.twin();
  }

  bool condition = false;

  bool is2Large = false;
  bool is2Long = false;
  bool is2Curved = false;
  bool is2Sharp = false;
  bool is2Fat = false;
  bool is2Skinny = false;
  bool isDelaunay = false;

  // const double targetdH0 = 0.5;

  // Conditions for splitting
  if (splitLarge) {
    is2Large = (isBoundary)
                   ? vpg.faceArea(he.face()) > (2 * targetFaceArea)
                   : (vpg.faceArea(he.face()) +
                      vpg.faceArea(he.twin().face())) > (4 * targetFaceArea);
    condition = condition || is2Large;
  }

  if (splitLong) {
    is2Long =
        (isBoundary)
            ? (vpg.edgeLength(e) > (2 * vpg.edgeLength(he.next().edge())))
            : (vpg.edgeLength(e) > (vpg.edgeLength(he.next().edge()) +
                                    vpg.edgeLength(he.twin().next().edge())));
    condition = is2Long || condition;
  }

  if (splitCurved) {
    is2Curved =
        vpg.edgeLength(e) > (2 * computeCurvatureThresholdLength(e, vpg));
    condition = is2Curved || condition;
  }

  if (splitSharp) {
    mem3dg_runtime_error("Split sharp is currently not supported");
    // is2Sharp = abs(H0[he.tipVertex()] - H0[he.tailVertex()]) > (2 *
    // targetdH0);
  }

  if (splitSkinnyDelaunay && !isBoundary) {
    isDelaunay = (vpg.cornerAngle(e.halfedge().next().next().corner()) +
                  vpg.cornerAngle(e.halfedge().twin().next().next().corner())) <
                 (constants::PI);
    is2Skinny =
        (vpg.cornerAngle(he.next().corner()) +
         vpg.cornerAngle(he.twin().next().corner())) < constants::PI / 3;
    condition = (isDelaunay && is2Skinny) || condition;
  }

  if (splitFat) {
    is2Fat = (isBoundary) ? (vpg.cornerAngle(he.next().next().corner()) >
                             constants::PI * 0.667)
                          : (vpg.cornerAngle(he.next().next().corner()) >
                                 constants::PI * 0.667 ||
                             vpg.cornerAngle(he.twin().next().next().corner()) >
                                 constants::PI * 0.667);
    condition = is2Fat || condition;
  }

  // bool flat = abs(vpg.edgeDihedralAngle(e)) < (constants::PI
  // / 36);
  // isSplit =
  //     is2Large || is2Curved || is2Fat || is2Sharp || (is2Skinny &&
  //     isDelaunay);
  return condition;
}

void MeshProcessor::MeshMutator::markVertices(
    gcs::VertexData<bool> &mutationMarker, const gcs::Vertex v,
    const size_t layer) {
  assert(("Max layer number is 2!", layer < 2));
  mutationMarker[v] = true;
  if (layer > 0) {
    for (gc::Vertex nv : v.adjacentVertices()) {
      mutationMarker[nv] = true;
      if (layer > 1) {
        for (gc::Vertex nnv : nv.adjacentVertices()) {
          mutationMarker[nnv] = true;
        }
      }
    }
  }
}

std::tuple<double, std::size_t> MeshProcessor::MeshMutator::neighborAreaSum(
    const gcs::Edge e, const gcs::VertexPositionGeometry &vpg) {
  double area = 0;
  std::size_t num_neighbor = -2;
  for (gcs::Vertex v : e.adjacentVertices()) {
    for (gcs::Face f : v.adjacentFaces()) {
      area += vpg.faceArea(f);
      num_neighbor++;
    }
  }
  return std::tie(area, num_neighbor);
}

double MeshProcessor::MeshMutator::computeCurvatureThresholdLength(
    const gcs::Edge e, const gcs::VertexPositionGeometry &vpg) {
  gcs::Halfedge he = e.halfedge();
  // curvature based remeshing:
  // https://www.irit.fr/recherches/VORTEX/publications/rendu-geometrie/EGshort2013_Dunyach_et_al.pdf
  auto tip = he.tipVertex();
  auto tail = he.tailVertex();
  double k1 = std::max(abs(vpg.vertexMaxPrincipalCurvature(tip)),
                       abs(vpg.vertexMinPrincipalCurvature(tip)));
  double k2 = std::max(abs(vpg.vertexMaxPrincipalCurvature(tail)),
                       abs(vpg.vertexMinPrincipalCurvature(tail)));

  return std::sqrt(6 * curvTol / std::max(k1, k2) - 3 * curvTol * curvTol);
}

} // namespace solver
} // namespace mem3dg
