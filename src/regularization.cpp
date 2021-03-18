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
        double &localReferenceArea = v.isBoundary()
                                         ? refVpg->faceAreas[base_he.face()]
                                         : meanTargetFaceArea;
        regularizationForce[v] +=
            -P.Ksl * localAreaGradient *
            (vpg->faceAreas[base_he.face()] - localReferenceArea);
      }

      // local edge regularization
      if (P.Kse != 0) {
        gc::Vector3 edgeGradient = -vecFromHalfedge(he, *vpg).normalize();
        double &localReferenceLength = v.isBoundary()
                                           ? refVpg->edgeLengths[he.edge()]
                                           : meanTargetEdgeLength;
        regularizationForce[v] +=
            -P.Kse * edgeGradient *
            (vpg->edgeLengths[he.edge()] - localReferenceLength);
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

  if (isFlipped) {
    mesh->compress();
  }

  return isFlipped;
}

template <typename E, typename T>
void averageData(gc::MeshData<E, T> &meshData, const E &element1, const E &element2,
                 const E &newElement) {
  meshData[newElement] = (meshData[element1] + meshData[element2]) / 2;
}
template <typename E, typename T>
T averageData(gc::MeshData<E, T> &meshData, const E &&element1, const E &&element2) {
	  return (meshData[element1] + meshData[element2]) / 2;
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
        gcs::Halfedge newhe = mesh->splitEdgeTriangular(he.edge());

        const auto &element1 = newhe.tipVertex(),
                   &element2 = newhe.next().next().twin().next().next().vertex(),
                   &newElement = newhe.vertex();
        averageData(vpg->inputVertexPositions, element1, element2, newElement);
        averageData(vel, element1, element2, newElement);
        if (O.isProtein)
          averageData(proteinDensity, element1, element2, newElement);

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

        auto meanPosition = averageData(vpg->inputVertexPositions,
                                        he.tipVertex(), he.tailVertex());
        auto meanVel = averageData(vel, he.tipVertex(), he.tailVertex());
        auto meanProteinDensity =
            averageData(proteinDensity, he.tipVertex(), he.tailVertex());

        auto newVertex = mesh->collapseEdgeTriangular(he.edge());

        if ((theVertex == he.tailVertex() || theVertex == he.tipVertex()) &&
            !O.isFloatVertex) // retrieve "the" vertex
          theVertex = newVertex;

        vpg->inputVertexPositions[newVertex] = meanPosition;
        vel[newVertex] = meanVel;
        if (O.isProtein)
          proteinDensity[newVertex] = meanProteinDensity;

        isGrown = true;
      }
    }
  }

  if (isGrown) {
    mesh->compress();
  }

  return isGrown;
}

} // namespace mem3dg
