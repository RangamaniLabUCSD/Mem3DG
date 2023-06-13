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
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"
#include "geometrycentral/utilities/vector3.h"
#include "mem3dg/constants.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/system.h"
#include <Eigen/Core>
#include <cmath>

#include <queue>

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

/**
 * @brief Check if an edge collapse will cause a foldover
 *
 * @param mesh
 * @param geometry
 * @param e
 * @return true
 * @return false
 */
bool System::ifFoldover(const gcs::Edge e, const double angle) const {
  gcs::Vertex v1 = e.firstVertex();
  gcs::Vertex v2 = e.secondVertex();

  const gc::Vector3 midpoint = (geometry.vpg->inputVertexPositions[v1] +
                                geometry.vpg->inputVertexPositions[v2]) /
                               2;

  // Check dihedral of triangle ABC and ABX to see if they majorly fold
  auto check = [&](gcs::Halfedge he0) {
    gcs::Halfedge heT = he0.twin();
    gc::Vector3 a = geometry.vpg->inputVertexPositions[heT.vertex()];
    gc::Vector3 b = geometry.vpg->inputVertexPositions[heT.next().vertex()];
    gc::Vector3 c =
        geometry.vpg->inputVertexPositions[heT.next().next().vertex()];

    gc::Vector3 n1 = gc::cross(b - a, c - a);
    gc::Vector3 n2 = gc::cross(b - midpoint, a - midpoint);
    return gc::angle(n1, n2) > angle;
  };

  // check link around v1
  for (gcs::Halfedge he : v1.outgoingHalfedges()) {
    auto candidate = he.next();
    if (candidate.vertex() != v2 && candidate.next().vertex() != v2) {
      if (check(candidate)) {
        return true;
      }
    }
  }

  // check link around v2
  for (gcs::Halfedge he : v2.outgoingHalfedges()) {
    auto candidate = he.next();
    if (candidate.vertex() != v1 && candidate.next().vertex() != v1) {
      if (check(candidate)) {
        return true;
      }
    }
  }

  return false;
}

void System::mutateMesh(size_t nRepetition) {
  for (size_t i = 0; i < nRepetition; ++i) {
    bool isGrown = false, isFlipped = false;
    mutationMarker.fill(false);

    // vertex shift for regularization
    if (meshProcessor.meshMutator.isShiftVertex) {
      vertexShift();
    }
    // linear edge flip for non-Delaunay triangles
    if (meshProcessor.meshMutator.isFlipEdge) {
      // isFlipped = edgeFlip();
      // isFlipped = edgeFlip() || isFlipped;
      // isFlipped = edgeFlip() || isFlipped;
      edgeFlip();
    }

    // split edge and collapse edge
    if (meshProcessor.meshMutator.isSplitEdge ||
        meshProcessor.meshMutator.isCollapseEdge) {
      isGrown = isGrown || processSplitCollapse();
    }

    // linear edge flip for non-Delaunay triangles
    if (meshProcessor.meshMutator.isFlipEdge) {
      // isFlipped = edgeFlip();
      // isFlipped = edgeFlip() || isFlipped;
      // isFlipped = edgeFlip() || isFlipped;
      edgeFlip();
    }

    if (meshProcessor.meshMutator.isSmoothenMesh) {
      smoothenMesh();
    }

    // globally update quantities
    if (isGrown || isFlipped) {
      globalUpdateAfterMutation();
    }
  }
}

void System::vertexShift() {
  for (gcs::Vertex v : geometry.mesh->vertices()) {
    // Only move if not significantly force masked
    if (gc::sum(forces.forceMask[v]) > 0.5) {
      if (v.isBoundary()) {
        gc::Vector3 barycenter{0.0, 0.0, 0.0};
        std::vector<gcs::Vertex> adjacent_boundary_vertices;
        adjacent_boundary_vertices.reserve(2);

        for (gcs::Vertex vAdj : v.adjacentVertices()) {
          if (vAdj.isBoundary()) {
            adjacent_boundary_vertices.push_back(vAdj);
          }
        }

        MEM3DG_SAFETY_ASSERT(
            adjacent_boundary_vertices.size() == 2,
            "Number of neighbor vertices on boundary is not 2!")

        const gcs::Vertex &v1 = adjacent_boundary_vertices[0];
        const gcs::Vertex &v2 = adjacent_boundary_vertices[1];

        barycenter = (geometry.vpg->inputVertexPositions[v1] +
                      geometry.vpg->inputVertexPositions[v2]) /
                     2;
        gc::Vector3 faceNormal =
            gc::cross(geometry.vpg->inputVertexPositions[v1] -
                          geometry.vpg->inputVertexPositions[v],
                      geometry.vpg->inputVertexPositions[v2] -
                          geometry.vpg->inputVertexPositions[v]);
        gc::Vector3 sideNormal =
            gc::cross(faceNormal, geometry.vpg->inputVertexPositions[v1] -
                                      geometry.vpg->inputVertexPositions[v2])
                .normalize();

        // project to move only in sidenormal direction
        geometry.vpg->inputVertexPositions[v] =
            barycenter -
            gc::dot(sideNormal,
                    barycenter - geometry.vpg->inputVertexPositions[v]) *
                sideNormal;
      } else {
        gc::Vector3 barycenter{0.0, 0.0, 0.0};

        int n_vAdj = 0;
        for (gcs::Vertex vAdj : v.adjacentVertices()) {
          barycenter += geometry.vpg->inputVertexPositions[vAdj];
          n_vAdj++;
        }
        barycenter /= n_vAdj;

        // geometry.vpg->inputVertexPositions[v] = barycenter;

        // for (gcs::Halfedge he : v.outgoingHalfedges()) {
        // gcs::Halfedge base_he = he.next();
        geometry.vpg->inputVertexPositions[v] =
            barycenter -
            gc::dot(geometry.vpg->vertexNormals[v],
                    barycenter - geometry.vpg->inputVertexPositions[v]) *
                geometry.vpg->vertexNormals[v];
        // }
      }
    }
  }
}

bool System::edgeFlip() {
  // Note in regularization, it is preferred to use immediate calculation rather
  // than cached one
  bool isFlipped = false;
  gcs::EdgeData<bool> isOrigEdge(*geometry.mesh, true);
  // flip edge if not delauney
  for (gcs::Edge e : geometry.mesh->edges()) {
    if (!isOrigEdge[e] || e.isBoundary()) {
      continue;
    }
    gcs::Halfedge he = e.halfedge();
    if (gc::sum(forces.forceMask[he.vertex()] +
                forces.forceMask[he.twin().vertex()]) < 0.5) {
      continue;
    }

    if (meshProcessor.meshMutator.checkFlipCondition(e, *geometry.vpg)) {
      bool success = geometry.mesh->flip(e);
      isOrigEdge[e] = false;
      isFlipped = true;
      meshProcessor.meshMutator.markVertices(mutationMarker, he.tailVertex());
      meshProcessor.meshMutator.markVertices(mutationMarker, he.tipVertex());
    }
  }

  if (isFlipped)
    geometry.mesh->compress();

  return isFlipped;
}

void System::edgeFlipQueued() {
  mem3dg_debug_message("Using queue based edge flip");

  // queue of edges to check if Delaunay
  std::queue<gc::Edge> toCheck;
  // true if edge is currently in toCheck
  gc::EdgeData<bool> inQueue(*geometry.mesh);
  // start with all edges
  for (gc::Edge e : geometry.mesh->edges()) {
    toCheck.push(e);
    inQueue[e] = true;
  }
  // Limiter for max flips
  int flipMax = 100 * geometry.mesh->nVertices();
  // Counter of flip attempts
  int flipCnt = 0;
  while (!toCheck.empty() && flipCnt < flipMax) {
    gc::Edge e = toCheck.front();
    toCheck.pop();
    inQueue[e] = false;
    // if not Delaunay, flip edge and enqueue the surrounding "diamond" edges
    // (if not already)
    if (meshProcessor.meshMutator.checkFlipCondition(e, *geometry.vpg)) {
      flipCnt++;
      gcs::Halfedge he = e.halfedge();

      // Dont flip if edge has masked force (enforce boundary condition)
      if (gc::sum(forces.forceMask[he.tailVertex()] +
                  forces.forceMask[he.tipVertex()]) < 0.5) {
        continue;
      }
      // Current triangle
      gcs::Halfedge he1 = he.next();
      gcs::Halfedge he2 = he1.next();
      // Twin triangle
      gcs::Halfedge he3 = he.twin().next();
      gcs::Halfedge he4 = he3.next();

      if (!inQueue[he1.edge()]) {
        toCheck.push(he1.edge());
        inQueue[he1.edge()] = true;
      }
      if (!inQueue[he2.edge()]) {
        toCheck.push(he2.edge());
        inQueue[he2.edge()] = true;
      }
      if (!inQueue[he3.edge()]) {
        toCheck.push(he3.edge());
        inQueue[he3.edge()] = true;
      }
      if (!inQueue[he4.edge()]) {
        toCheck.push(he4.edge());
        inQueue[he4.edge()] = true;
      }
      geometry.mesh->flip(e);

      meshProcessor.meshMutator.markVertices(mutationMarker, he.tailVertex());
      meshProcessor.meshMutator.markVertices(mutationMarker, he.tipVertex());
    }
  }
}

bool System::processSplitCollapse() {
  // Note in regularization, it is preferred to use immediate calculation rather
  // than cached one
  bool isGrown = false;
  gcs::EdgeData<bool> isOrigEdge(*geometry.mesh, true);
  // gcs::VertexData<bool> isOrigVertex(*geometry.mesh, true);

  // expand the mesh when area is too large
  for (gcs::Edge e : geometry.mesh->edges()) {

    // don't keep processing new edges
    if (!isOrigEdge[e])
      continue;

    // alias the halfedge
    gcs::Halfedge he = e.halfedge();

    // gather both vertices and their properties
    gcs::Vertex vertex1 = he.tipVertex(), vertex2 = he.tailVertex();
    gc::Vector3 vertex1ForceMask = forces.forceMask[vertex1];
    gc::Vector3 vertex2ForceMask = forces.forceMask[vertex2];
    // don't keep processing static vertices
    if (gc::sum(vertex1ForceMask + vertex2ForceMask) < 0.5)
      continue;

    // Splitting
    if (meshProcessor.meshMutator.checkSplitCondition(e, *geometry.vpg)) {
      auto newVertex = splitEdge(e);

      // isOrigVertex[newVertex] = false;
      for (gcs::Edge adjEdge : newVertex.adjacentEdges()) {
        isOrigEdge[adjEdge] = false;
      }

      meshProcessor.meshMutator.markVertices(mutationMarker, newVertex);
      // mutationMarker[newVertex] = true;

      isGrown = true;
    } else if (meshProcessor.meshMutator.checkCollapseCondition(
                   e, *geometry.vpg) &&
               !ifFoldover(e)) { // Collapsing
      auto newVertex = collapseEdge(e);
      if (newVertex.getIndex() != gc::INVALID_IND) {
        // isOrigVertex[newVertex] = false;
        for (gcs::Edge e : newVertex.adjacentEdges()) {
          isOrigEdge[e] = false;
        }
        meshProcessor.meshMutator.markVertices(mutationMarker, newVertex);
        isGrown = true;
      }
    }
  } // end for edge
  if (isGrown) {
    geometry.mesh->compress();
  }
  return isGrown;
}

bool System::processSplitCollapseQueued() {
  bool didSplitOrCollapse = false;
  // queues of edges to CHECK to change
  std::vector<gcs::Edge> toSplit;
  std::vector<gcs::Edge> toCollapse;

  for (gc::Edge e : geometry.mesh->edges()) {
    toSplit.push_back(e);
  }

  // actually splitting
  while (!toSplit.empty()) {
    gcs::Edge e = toSplit.back();
    toSplit.pop_back();

    // gcs::Halfedge he = e.halfedge();

    gcs::Vertex vertex1 = e.firstVertex(), vertex2 = e.secondVertex();
    gc::Vector3 vertex1ForceMask = forces.forceMask[vertex1];
    gc::Vector3 vertex2ForceMask = forces.forceMask[vertex2];

    if (meshProcessor.meshMutator.checkSplitCondition(e, *geometry.vpg) &&
        gc::sum(vertex1ForceMask + vertex2ForceMask) > 0.5) {

      auto newVertex = splitEdge(e);

      didSplitOrCollapse = true;
      meshProcessor.meshMutator.markVertices(mutationMarker, newVertex);
    } else {
      toCollapse.push_back(e);
    }
  }
  // actually collapsing
  while (!toCollapse.empty()) {
    gcs::Edge e = toCollapse.back();
    toCollapse.pop_back();
    if (e.halfedge().next().getIndex() !=
        gc::INVALID_IND) { // make sure it exists
      gcs::Halfedge he = e.halfedge();
      gcs::Vertex vertex1 = he.tipVertex(), vertex2 = he.tailVertex();
      gc::Vector3 vertex1ForceMask = forces.forceMask[vertex1];
      gc::Vector3 vertex2ForceMask = forces.forceMask[vertex2];

      if (meshProcessor.meshMutator.checkCollapseCondition(e, *geometry.vpg) &&
          gc::sum(vertex1ForceMask + vertex2ForceMask) > 0.5 &&
          !ifFoldover(e)) {

        auto newVertex = collapseEdge(e);
        if (newVertex.getIndex() != gc::INVALID_IND) {
          meshProcessor.meshMutator.markVertices(mutationMarker, newVertex);
          didSplitOrCollapse = true;
        }
      }
    }
  }
  if (didSplitOrCollapse)
    geometry.mesh->compress();
  return didSplitOrCollapse;
}

Eigen::Matrix<bool, Eigen::Dynamic, 1>
System::smoothenMesh(double initStep, double target, size_t maxIteration) {
  // require nonzero bending rigidity in parameters
  if (Kb.raw().sum() == 0) {
    mem3dg_runtime_error(
        "Bending rigidity has to be nonzero to smoothen mesh!");
  }
  // initialize variables
  double stepSize = initStep;
  double pastGradNorm = 1e10;
  size_t num_iter = 0;
  // compute bending forces
  geometry.vpg->refreshQuantities();
  computeGeometricForces();
  EigenVectorX3dr pastForceVec(toMatrix(forces.spontaneousCurvatureForceVec));

  // initialize smoothingMask: we only smooth if there are outliers in the
  // spontaneous curvature force
  Eigen::Matrix<bool, Eigen::Dynamic, 1> smoothingMask =
      outlierMask(forces.spontaneousCurvatureForce.raw(), 0.5);
  isSmooth =
      (std::count(smoothingMask.begin(), smoothingMask.end(), true) == 0);
  // initialize gradient and compute exit tolerance
  double gradNorm = toMatrix(forces.spontaneousCurvatureForceVec).norm();
  double tol = gradNorm * target;

  while (gradNorm > tol && !isSmooth) {
    if (stepSize < 1e-8 * initStep) {
      mem3dg_runtime_warning("smoothing operation diverges: ", stepSize,
                             " < 1e-8.");
    }
    if (num_iter == maxIteration) {
      mem3dg_runtime_warning("smoothing operation exceeds max iteration!");
      break;
    }

    // compute spontaneous curvature force if smoothingMask is true
    geometry.vpg->refreshQuantities();
    forces.spontaneousCurvatureForceVec.fill({0, 0, 0});
    forces.spontaneousCurvatureForce.raw().setZero();
    for (std::size_t i = 0; i < geometry.mesh->nVertices(); ++i) {
      if (smoothingMask[i]) {
        computeGeometricForces(i);
      }
    }
    // compute norm of the spontaneous curvature force
    gradNorm = toMatrix(forces.spontaneousCurvatureForceVec).norm();
    // Previous step has led to an increase in force. Recover the position and
    // cut the step size in half
    if (gradNorm > pastGradNorm) {
      toMatrix(geometry.vpg->inputVertexPositions) -= pastForceVec * stepSize;
      stepSize /= 2;
      continue;
    }
    // smoothing step
    geometry.vpg->inputVertexPositions +=
        forces.spontaneousCurvatureForceVec * stepSize;
    pastGradNorm = gradNorm;
    pastForceVec = toMatrix(forces.spontaneousCurvatureForceVec);
    num_iter++;
  };

  return smoothingMask;
}

void System::localSmoothing(const gcs::Vertex &v, std::size_t num,
                            double stepSize) {
  std::size_t count = 0;
  while (count < num) {
    gc::Vector3 vertexNormal{0, 0, 0};
    for (gcs::Corner c : v.adjacentCorners()) {
      vertexNormal +=
          geometry.vpg->cornerAngle(c) * geometry.vpg->faceNormal(c.face());
    }
    vertexNormal.normalize();
    double localLapH = 0;
    double H_center =
        geometry.vpg->vertexMeanCurvature(v) / geometry.vpg->vertexDualArea(v);
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      localLapH +=
          geometry.vpg->edgeCotanWeight(he.edge()) *
          (H_center - geometry.vpg->vertexMeanCurvature(he.tipVertex()) /
                          geometry.vpg->vertexDualArea(he.tipVertex()));
    }
    geometry.vpg->inputVertexPositions[v] -=
        stepSize * geometry.vpg->vertexDualArea(v) * localLapH * vertexNormal;
    count++;
  }
}

void System::localSmoothing(const gcs::Halfedge &he, std::size_t num,
                            double stepSize) {
  std::size_t count = 0;
  while (count < num) {
    gc::Vector3 vertexNormal1{0, 0, 0};
    gc::Vector3 vertexNormal2{0, 0, 0};
    double localLapH1 = 0;
    double localLapH2 = 0;

    gcs::Vertex v = he.tailVertex();
    for (gcs::Corner c : v.adjacentCorners()) {
      vertexNormal1 +=
          geometry.vpg->cornerAngle(c) * geometry.vpg->faceNormal(c.face());
    }
    vertexNormal1.normalize();
    double H_center =
        geometry.vpg->vertexMeanCurvature(v) / geometry.vpg->vertexDualArea(v);
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      localLapH1 +=
          geometry.vpg->edgeCotanWeight(he.edge()) *
          (H_center - geometry.vpg->vertexMeanCurvature(he.tipVertex()) /
                          geometry.vpg->vertexDualArea(he.tipVertex()));
    }

    v = he.tipVertex();
    for (gcs::Corner c : v.adjacentCorners()) {
      vertexNormal2 +=
          geometry.vpg->cornerAngle(c) * geometry.vpg->faceNormal(c.face());
    }
    vertexNormal2.normalize();
    H_center =
        geometry.vpg->vertexMeanCurvature(v) / geometry.vpg->vertexDualArea(v);
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      localLapH2 +=
          geometry.vpg->edgeCotanWeight(he.edge()) *
          (H_center - geometry.vpg->vertexMeanCurvature(he.tipVertex()) /
                          geometry.vpg->vertexDualArea(he.tipVertex()));
    }

    geometry.vpg->inputVertexPositions[he.tailVertex()] -=
        stepSize * geometry.vpg->vertexDualArea(he.tailVertex()) * localLapH1 *
        vertexNormal1;
    geometry.vpg->inputVertexPositions[he.tipVertex()] -=
        stepSize * geometry.vpg->vertexDualArea(he.tipVertex()) * localLapH2 *
        vertexNormal2;
    count++;
  }
}

void System::globalUpdateAfterMutation() {
  // update the velocity
  velocity = forces.maskForce(velocity); // important: velocity interpolation
                                         // contaminate the zero velocity
  if (computeKineticEnergy() != 0) {
    double oldKE = energy.kineticEnergy;
    velocity *= pow(oldKE / computeKineticEnergy(), 0.5);
  }

  // Update mask when topology changes (likely not necessary, just for safety)
  if (geometry.mesh->hasBoundary()) {
    forces.forceMask.fill({1, 1, 1});
    boundaryForceMask(*geometry.mesh, forces.forceMask,
                      parameters.boundary.shapeBoundaryCondition);
    forces.proteinMask.fill(1);
    boundaryProteinMask(*geometry.mesh, forces.proteinMask,
                        parameters.boundary.proteinBoundaryCondition);
  }

  // if (geometry.notableVertex.raw().cast<int>().sum() != 1)
  //   mem3dg_runtime_error("there are more than one true in center!");
}

} // namespace solver
} // namespace mem3dg
