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
#include <cmath>

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void MeshProcessor::summarizeStatus() {
  meshRegularizer.summarizeStatus();
  meshMutator.summarizeStatus();

  isMeshRegularize = (meshRegularizer.Kst != 0) || (meshRegularizer.Ksl != 0) ||
                     (meshRegularizer.Kse != 0);
  isMeshMutate = meshMutator.isChangeTopology || meshMutator.shiftVertex;

  if (!meshRegularizer.ifHasRefMesh && isMeshRegularize) {
    mem3dg_runtime_error(
        "To apply mesh regularization, reference mesh has to be provided!");
  }
  if (meshMutator.isChangeTopology && isMeshRegularize) {
    mem3dg_runtime_error("For topology changing simulation, mesh "
                         "regularization cannot be applied!");
  }
};

void MeshProcessor::MeshRegularizer::summarizeStatus() {
  ifHasRefMesh = (nEdge != 0) || (nVertex != 0) || (nFace != 0);
};

void MeshProcessor::MeshRegularizer::readReferenceData(std::string refMesh,
                                                       std::size_t nSub) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> referenceMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> referenceVpg;

  // Load input reference mesh and geometry
  std::tie(referenceMesh, referenceVpg) = gcs::readManifoldSurfaceMesh(refMesh);

  // Subdivide the mesh and geometry objects
  if (nSub > 0) {
    // mem3dg::subdivide(ptrRefMesh, ptrRefVpg, nSub);
    mem3dg::loopSubdivide(referenceMesh, referenceVpg, nSub);
  }

  nEdge = referenceMesh->nEdges();
  nVertex = referenceMesh->nVertices();
  nFace = referenceMesh->nFaces();

  referenceVpg->requireEdgeLengths();
  referenceVpg->requireFaceAreas();
  refFaceAreas = referenceVpg->faceAreas.raw();
  refEdgeLengths = referenceVpg->edgeLengths.raw();
  refLcrs.resize(referenceMesh->nEdges(), 1);
  for (std::size_t i = 0; i < referenceMesh->nEdges(); ++i) {
    gcs::Edge e{referenceMesh->edge(i)};
    refLcrs[i] = computeLengthCrossRatio(*referenceVpg, e);
  }

  meanTargetFaceArea = refFaceAreas.sum() / referenceMesh->nFaces();
  meanTargetEdgeLength = refEdgeLengths.sum() / referenceMesh->nEdges();
}

void MeshProcessor::MeshRegularizer::readReferenceData(
    Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &topologyMatrix,
    Eigen::Matrix<double, Eigen::Dynamic, 3> &refVertexMatrix,
    std::size_t nSub) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> referenceMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> referenceVpg;

  // Load input reference mesh and geometry
  std::tie(referenceMesh, referenceVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(refVertexMatrix, topologyMatrix);
  std::cout << "Loaded reference mesh " << std::endl;

  // Subdivide the mesh and geometry objects
  if (nSub > 0) {
    // mem3dg::subdivide(ptrRefMesh, ptrRefVpg, nSub);
    mem3dg::loopSubdivide(referenceMesh, referenceVpg, nSub);
  }

  nEdge = referenceMesh->nEdges();
  nVertex = referenceMesh->nVertices();
  nFace = referenceMesh->nFaces();

  referenceVpg->requireEdgeLengths();
  referenceVpg->requireFaceAreas();
  refFaceAreas = referenceVpg->faceAreas.raw();
  refEdgeLengths = referenceVpg->edgeLengths.raw();
  refLcrs.resize(referenceMesh->nEdges(), 1);
  for (std::size_t i = 0; i < referenceMesh->nEdges(); ++i) {
    gcs::Edge e{referenceMesh->edge(i)};
    refLcrs[i] = computeLengthCrossRatio(*referenceVpg, e);
  }

  meanTargetFaceArea = refFaceAreas.sum() / referenceMesh->nFaces();
  meanTargetEdgeLength = refEdgeLengths.sum() / referenceMesh->nEdges();
}

double MeshProcessor::MeshRegularizer::computeLengthCrossRatio(
    gcs::VertexPositionGeometry &vpg, gcs::Edge &e) const {

  gcs::Edge lj = e.halfedge().next().edge();
  gcs::Edge ki = e.halfedge().twin().next().edge();
  gcs::Edge il = e.halfedge().next().next().edge();
  gcs::Edge jk = e.halfedge().twin().next().next().edge();
  return vpg.edgeLengths[il] * vpg.edgeLengths[jk] / vpg.edgeLengths[ki] /
         vpg.edgeLengths[lj];
}

bool MeshProcessor::MeshMutator::ifFlip(
    const gcs::Edge e, const gcs::VertexPositionGeometry &vpg) {
  gcs::Halfedge he = e.halfedge();
  bool condition = false;

  if (flipNonDelaunay && !e.isBoundary()) {
    bool nonDelaunay =
        (vpg.cornerAngle(he.next().next().corner()) +
         vpg.cornerAngle(he.twin().next().next().corner())) > (constants::PI);
    if (flipNonDelaunayRequireFlat) {
      bool flat = abs(vpg.edgeDihedralAngle(he.edge())) < (constants::PI / 36);
      nonDelaunay = nonDelaunay && flat;
    }
    condition = condition || nonDelaunay;
  }

  return condition;
}

void MeshProcessor::MeshMutator::summarizeStatus() {
  isEdgeFlip = (flipNonDelaunay || flipNonDelaunayRequireFlat);
  isSplitEdge = (splitCurved || splitLarge || splitLong || splitSharp ||
                 splitSkinnyDelaunay);
  isCollapseEdge = (collapseSkinny || collapseSmall || collapseSmallNeedFlat);
  isChangeTopology = isEdgeFlip || isSplitEdge || isCollapseEdge;
};

bool MeshProcessor::MeshMutator::ifCollapse(
    const gc::Edge e, const gcs::VertexPositionGeometry &vpg) {
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
    is2Skinny =
        (isBoundary)
            ? (vpg.cornerAngle(he.next().next().corner()) < constants::PI / 6)
            : (vpg.cornerAngle(he.next().next().corner()) +
               vpg.cornerAngle(he.twin().next().next().corner())) <
                      constants::PI / 3 ||
                  (vpg.cornerAngle(he.next().corner()) +
                   vpg.cornerAngle(he.twin().corner())) >
                      (constants::PI * 1.333) ||
                  (vpg.cornerAngle(he.corner()) +
                   vpg.cornerAngle(he.twin().next().corner())) >
                      (constants::PI * 1.333);
    condition = is2Skinny;
  }

  if (collapseSmall) {
    double areaSum;
    std::size_t num_neighbor;
    std::tie(areaSum, num_neighbor) = neighborAreaSum(e, vpg);

    size_t gap = isBoundary ? 1 : 2;
    is2Small =
        (areaSum - vpg.faceArea(he.face()) - vpg.faceArea(he.twin().face())) <
        (num_neighbor - gap) * targetFaceArea;
    if (collapseSmallNeedFlat) {
      isFlat =
          vpg.edgeLength(e) < (0.667 * computeCurvatureThresholdLength(e, vpg));
      is2Small = is2Small && isFlat;
    }
    // isSmooth =
    //     abs(H0[he.tipVertex()] - H0[he.tailVertex()]) < (0.667 *
    //     targetdH0);
    condition = condition || is2Small;
  }

  // isCollapse = is2Skinny; //|| (is2Small && isFlat && isSmooth);

  return condition;
}

bool MeshProcessor::MeshMutator::ifSplit(
    const gcs::Edge e, const gcs::VertexPositionGeometry &vpg) {
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
  if (layer > 2)
    mem3dg_runtime_error("max layer number is 2!");
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
  double k1 = (abs(vpg.vertexMaxPrincipalCurvature(he.tipVertex())) >
               abs(vpg.vertexMinPrincipalCurvature(he.tipVertex())))
                  ? abs(vpg.vertexMaxPrincipalCurvature(he.tipVertex()))
                  : abs(vpg.vertexMinPrincipalCurvature(he.tipVertex()));
  double k2 = (abs(vpg.vertexMaxPrincipalCurvature(he.tailVertex())) >
               abs(vpg.vertexMinPrincipalCurvature(he.tailVertex())))
                  ? abs(vpg.vertexMaxPrincipalCurvature(he.tailVertex()))
                  : abs(vpg.vertexMinPrincipalCurvature(he.tailVertex()));
  return std::sqrt(6 * curvTol / ((k1 > k2) ? k1 : k2) - 3 * curvTol * curvTol);
}

} // namespace solver
} // namespace mem3dg
