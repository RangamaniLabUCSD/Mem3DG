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
#include "mem3dg/solver/constants.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"
#include <Eigen/Core>

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

gcs::EdgeData<double>
System::computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg) const {
  gcs::EdgeData<double> LCR(*mesh);
  for (gcs::Edge e : mesh->edges()) {
    gcs::Edge lj = e.halfedge().next().edge();
    gcs::Edge ki = e.halfedge().twin().next().edge();
    gcs::Edge il = e.halfedge().next().next().edge();
    gcs::Edge jk = e.halfedge().twin().next().next().edge();
    // clr[e] = edgeLength[il] * edgeLength[jk] / edgeLength[ki] /
    //          edgeLength[lj];
    LCR[e] = vpg.edgeLengths[il] * vpg.edgeLengths[jk] / vpg.edgeLengths[ki] /
             vpg.edgeLengths[lj];
  }
  return LCR;
}

double System::computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                                       gcs::Edge &e) const {

  gcs::Edge lj = e.halfedge().next().edge();
  gcs::Edge ki = e.halfedge().twin().next().edge();
  gcs::Edge il = e.halfedge().next().next().edge();
  gcs::Edge jk = e.halfedge().twin().next().next().edge();
  return vpg.edgeLengths[il] * vpg.edgeLengths[jk] / vpg.edgeLengths[ki] /
         vpg.edgeLengths[lj];
}

double System::computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                                       gcs::Edge &&e) const {

  gcs::Edge lj = e.halfedge().next().edge();
  gcs::Edge ki = e.halfedge().twin().next().edge();
  gcs::Edge il = e.halfedge().next().next().edge();
  gcs::Edge jk = e.halfedge().twin().next().next().edge();
  return vpg.edgeLengths[il] * vpg.edgeLengths[jk] / vpg.edgeLengths[ki] /
         vpg.edgeLengths[lj];
}

void System::computeRegularizationForce() {
  for (gcs::Vertex v : mesh->vertices()) {
    for (gcs::Halfedge he : v.outgoingHalfedges()) {

      // Conformal regularization
      if (P.Kst != 0 && !he.edge().isBoundary()) {
        gcs::Halfedge jl = he.next();
        gcs::Halfedge li = jl.next();
        gcs::Halfedge ik = he.twin().next();
        gcs::Halfedge kj = ik.next();

        gc::Vector3 grad_li = vecFromHalfedge(li, *vpg).normalize();
        gc::Vector3 grad_ik = vecFromHalfedge(ik.twin(), *vpg).normalize();
        regularizationForce[v] +=
            -P.Kst *
            (computeLengthCrossRatio(*vpg, he.edge()) - targetLcr[he.edge()]) /
            targetLcr[he.edge()] *
            (vpg->edgeLengths[kj.edge()] / vpg->edgeLengths[jl.edge()]) *
            (grad_li * vpg->edgeLengths[ik.edge()] -
             grad_ik * vpg->edgeLengths[li.edge()]) /
            vpg->edgeLengths[ik.edge()] / vpg->edgeLengths[ik.edge()];
      }

      // Local area regularization
      if (P.Ksl != 0 && he.isInterior()) {
        gcs::Halfedge base_he = he.next();
        gc::Vector3 base_vec = vecFromHalfedge(base_he, *vpg);
        gc::Vector3 localAreaGradient =
            -gc::cross(base_vec, vpg->faceNormals[he.face()]);
        // ->faceArea() is used over ->faceAreas[] when initConst() is not on
        // refVpg
        regularizationForce[v] +=
            -P.Ksl * localAreaGradient *
            (vpg->faceAreas[base_he.face()] - v.isBoundary()
                 ? refVpg->faceArea(base_he.face())
                 : meanTargetFaceArea);
      }

      // local edge regularization
      if (P.Kse != 0) {
        gc::Vector3 edgeGradient = -vecFromHalfedge(he, *vpg).normalize();
        regularizationForce[v] += -P.Kse * edgeGradient *
                                  (vpg->edgeLengths[he.edge()] - v.isBoundary()
                                       ? refVpg->edgeLength(he.edge())
                                       : meanTargetEdgeLength);
      }
    }
  }

  // post processing regularization force
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg->vertexNormals);
  auto regularizationForce_e = gc::EigenMap<double, 3>(regularizationForce);

  // remove the normal component
  regularizationForce_e -= rowwiseScaling(
      rowwiseDotProduct(regularizationForce_e, vertexAngleNormal_e),
      vertexAngleNormal_e);

  // moving boundary
  for (gcs::Vertex v : mesh->vertices()) {
    if (!mask[v]) {
      regularizationForce[v].z = 0;
    }
  }

  // / Patch regularization
  // the cubic penalty is for regularizing the mesh,
  // need better physical interpretation or alternative method
  // if (P.Kse != 0) {
  //   double strain =
  //       (vpg->edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
  //       targetEdgeLengths[he.edge()];
  //   regularizationForce[v] +=
  //      -P.Kse * edgeGradient * strain * strain * strain;
  // }

  // // remove the masked components
  // regularizationForce_e =
  //     rowwiseScaling(mask.raw().cast<double>(), regularizationForce_e);
  // // moving boundary
  // for (gcs::Vertex v : mesh->vertices()) {
  //   if (!mask[v]) {
  //     regularizationForce[v].z = 0;
  //     // boundary tension, mostly likely not necessary
  //     if (v.isBoundary()) {
  //       double boundaryEdgeLength = 0;
  //       for (gcs::Edge e : v.adjacentEdges()) {
  //         if (e.isBoundary()) {
  //           boundaryEdgeLength += vpg->edgeLength(e);
  //         }
  //       }
  //       boundaryEdgeLength /= 2;
  //       double scaling;
  //       if (regularizationForce[v].norm() > 1e-15) {
  //         scaling = 1 - abs(P.Ksg * boundaryEdgeLength /
  //                           (regularizationForce[v].norm()));
  //       }
  //       regularizationForce[v] *= scaling;
  //     }
  //   }
  // }
}

void System::vertexShift() {
  for (gcs::Vertex v : mesh->vertices()) {
    if (mask[v]) {
      gc::Vector3 baryCenter{0.0, 0.0, 0.0};
      double n_vAdj = 0.0;
      for (gcs::Vertex vAdj : v.adjacentVertices()) {
        baryCenter += vpg->inputVertexPositions[vAdj];
        n_vAdj += 1.0;
      }
      baryCenter /= n_vAdj;
      for (gcs::Halfedge he : v.outgoingHalfedges()) {
        gcs::Halfedge base_he = he.next();
        vpg->inputVertexPositions[v] =
            baryCenter - gc::dot(vpg->vertexNormals[v],
                                 baryCenter - vpg->inputVertexPositions[v]) *
                             vpg->vertexNormals[v];
      }
    }
  }
}

bool System::edgeFlip() {
  bool isFlipped = false;

  // flip edge if not delauney
  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    if (mask[he.vertex()] && mask[he.twin().vertex()]) {
      if ((vpg->cornerAngle(he.next().next().corner()) +
           vpg->cornerAngle(he.twin().next().next().corner())) >
          constants::PI) {
        // isFlip[he.edge()] = true;
        auto success = mesh->flip(he.edge());
        isFlipped = true;
      }
    }
  }

  if (isFlipped)
    mesh->compress();

  return isFlipped;
}

bool System::growMesh() {
  bool isGrown = false;
  // expand the mesh when area is too large
  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    if (mask[he.vertex()] && mask[he.twin().vertex()]) {
      // if(!he.edge().isBoundary()){
      if ((vpg->faceArea(he.face()) + vpg->faceArea(he.twin().face())) >
          (4 * meanTargetFaceArea)) {
        // || abs(H0[he.tailVertex()] - H0[he.tipVertex()]) > 0.2) {
        // split the edge
        const auto &vertex1 = he.tipVertex(), &vertex2 = he.tailVertex();
        gcs::Halfedge &newhe = mesh->splitEdgeTriangular(he.edge());
        const auto &newVertex = newhe.vertex();

        // update quantities
        localUpdateAfterMutation(vertex1, vertex2, newVertex);

        isGrown = true;
      }
    }
  }

  // shrink the mesh is vertex is too close
  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    if (mask[he.vertex()] && mask[he.twin().vertex()]) {
      if ((vpg->cornerAngle(he.next().next().corner()) +
           vpg->cornerAngle(he.twin().next().next().corner())) <
          constants::PI / 3) {
        // collapse the edge
        const auto &vertex1 = he.tipVertex(), &vertex2 = he.tailVertex();
        auto newVertex = mesh->collapseEdgeTriangular(he.edge());

        // update quantities
        localUpdateAfterMutation(vertex1, vertex2, newVertex);

        // update "the" vertex
        // if (P.pt.size() == 1 &&
        //     (thePoint.vertex == vertex1 || thePoint.vertex == vertex2) &&
        //     !O.isFloatVertex)
        //   thePoint.vertex = newVertex;

        isGrown = true;
      }
    }
  }

  if (isGrown)
    mesh->compress();
  return isGrown;
}

void System::processMesh() {

  bool isTopologyChanged = false;

  // vertex shift for regularization
  if (O.isVertexShift) {
    vertexShift();
  }

  // split edge and collapse edge
  if (O.isGrowMesh) {
    isTopologyChanged = growMesh();
  }

  // linear edge flip for non-Delauney triangles
  if (O.isEdgeFlip) {
    isTopologyChanged = edgeFlip();
  }

  // regularization
  if ((P.Kse != 0) || (P.Ksl != 0) || (P.Kst != 0)) {
    computeRegularizationForce();
    vpg->inputVertexPositions.raw() += regularizationForce.raw();
  }

  // globally update quantities
  globalUpdateAfterMeshProcessing(isTopologyChanged);
}

template <typename T>
void System::localUpdateAfterMutation(const T &element1, const T &element2,
                                      const T &newElement) {
  averageData(vpg->inputVertexPositions, element1, element2, newElement);
  averageData(vel, element1, element2, newElement);
  if (O.isProtein)
    averageData(proteinDensity, element1, element2, newElement);
}

void System::globalUpdateAfterMeshProcessing(bool &isTopologyChanged) {

  // recompute "the vertex" after topological changes
  if (O.isFloatVertex) {
    findTheVertex(*vpg);
  }

  // initialize/update local spontaneous curvature and bending rigidity
  if (O.isLocalCurvature || P.Kf != 0) {
    vpg->refreshQuantities();

    if (O.isGrowMesh || O.isEdgeFlip) {
      gcs::HeatMethodDistanceSolver heatSolverLocal(*vpg);
      geodesicDistanceFromPtInd = heatSolverLocal.computeDistance(thePoint);
    } else {
      geodesicDistanceFromPtInd = heatSolver.computeDistance(thePoint);
    }

    tanhDistribution(*vpg, H0.raw(), geodesicDistanceFromPtInd.raw(),
                     P.sharpness, P.r_H0);
    tanhDistribution(*vpg, Kb.raw(), geodesicDistanceFromPtInd.raw(),
                     P.sharpness, P.r_H0);
    H0.raw() *= P.H0;
    Kb.raw() *= P.Kbc - P.Kb;
    Kb.raw().array() += P.Kb;
    // ellipticDistribution(*vpg, H0.raw(), geodesicDistanceFromPtInd.raw(),
    //                      P.r_H0);
    // ellipticDistribution(*vpg, Kb.raw(), geodesicDistanceFromPtInd.raw(),
    //                      P.r_H0);
  }

  if (isTopologyChanged) {
    // Update the distribution matrix when topology changes
    if (P.eta != 0) {
      D = vpg->d0.transpose();
      for (int k = 0; k < D.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it) {
          it.valueRef() = 0.5;
        }
      }
    }

    // Update mask when topology changes
    if (O.isOpenMesh) {
      boundaryMask(*mesh, mask.raw());
    }

    // Update spontaneous curvature and bending rigidity when topology changes
    if (!O.isLocalCurvature) {
      H0.raw().setConstant(mesh->nVertices(), 1, P.H0);
      Kb.raw().setConstant(mesh->nVertices(), 1, P.Kb);
    }
  }
}

} // namespace mem3dg
