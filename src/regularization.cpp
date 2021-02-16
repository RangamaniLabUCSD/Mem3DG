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
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"
#include <Eigen/Core>

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

gcs::EdgeData<double> System::getLengthCrossRatio(gcs::VertexPositionGeometry &vpg) const{
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

void System::getRegularizationForce() {
  gcs::EdgeData<double> lcr(*mesh);
  lcr = getLengthCrossRatio(*vpg);

  for (gcs::Vertex v : mesh->vertices()) {
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gcs::Halfedge base_he = he.next();

      // Stretching forces
      gc::Vector3 edgeGradient = -vecFromHalfedge(he, *vpg).normalize();
      gc::Vector3 base_vec = vecFromHalfedge(base_he, *vpg);
      gc::Vector3 localAreaGradient =
          -gc::cross(base_vec, vpg->faceNormals[he.face()]);
      assert((gc::dot(localAreaGradient, vecFromHalfedge(he, *vpg))) < 0);

      // Conformal regularization
      if (P.Kst != 0) {
        gcs::Halfedge jl = he.next();
        gcs::Halfedge li = jl.next();
        gcs::Halfedge ik = he.twin().next();
        gcs::Halfedge kj = ik.next();

        gc::Vector3 grad_li = vecFromHalfedge(li, *vpg).normalize();
        gc::Vector3 grad_ik = vecFromHalfedge(ik.twin(), *vpg).normalize();
        regularizationForce[v] +=
            -P.Kst * (lcr[he.edge()] - targetLcr[he.edge()]) /
            targetLcr[he.edge()] *
            (vpg->edgeLengths[kj.edge()] / vpg->edgeLengths[jl.edge()]) *
            (grad_li * vpg->edgeLengths[ik.edge()] -
             grad_ik * vpg->edgeLengths[li.edge()]) /
            vpg->edgeLengths[ik.edge()] / vpg->edgeLengths[ik.edge()];
        // regularizationForce[v] += -P.Kst * localAreaGradient;
      }

      if (P.Ksl != 0) {
        regularizationForce[v] +=
            -P.Ksl * localAreaGradient *
            (vpg->faceAreas[base_he.face()] - targetFaceAreas[base_he.face()]) /
            targetFaceAreas[base_he.face()];
      }

      if (P.Kse != 0) {
        regularizationForce[v] +=
            -P.Kse * edgeGradient *
            (vpg->edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
            targetEdgeLengths[he.edge()];
      }

      /// Patch regularization
      // // the cubic penalty is for regularizing the mesh,
      // // need better physical interpretation or alternative method
      // if (P.Kse != 0) {
      //   double strain =
      //       (vpg->edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
      //       targetEdgeLengths[he.edge()];
      //   regularizationForce[v] +=
      //       -P.Kse * edgeGradient * strain * strain * strain;
      // }
    }
  }
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

void System::edgeFlip() {
}

void System::growMesh() {

}

} // namespace mem3dg