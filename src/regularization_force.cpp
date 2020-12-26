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

#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/IterativeLinearSolvers>

#include "mem3dg/solver/force.h"
#include "mem3dg/solver/util.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::getRegularizationForce() {
  gcs::EdgeData<double> lcr(mesh);
  getCrossLengthRatio(mesh, vpg, lcr);
  
  for (gcs::Vertex v : mesh.vertices()) {
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gcs::Halfedge base_he = he.next();

      // Stretching forces
      gc::Vector3 edgeGradient = -vecFromHalfedge(he, vpg).normalize();
      gc::Vector3 base_vec = vecFromHalfedge(base_he, vpg);
      gc::Vector3 localAreaGradient =
          -gc::cross(base_vec, vpg.faceNormals[he.face()]);
      assert((gc::dot(localAreaGradient, vecFromHalfedge(he, vpg))) < 0);

      // Conformal regularization
      if (P.Kst != 0) {
        gcs::Halfedge jl = he.next();
        gcs::Halfedge li = jl.next();
        gcs::Halfedge ik = he.twin().next();
        gcs::Halfedge kj = ik.next();

        gc::Vector3 grad_li = vecFromHalfedge(li, vpg).normalize();
        gc::Vector3 grad_ik = vecFromHalfedge(ik.twin(), vpg).normalize();
        regularizationForce[v] +=
            -P.Kst * (lcr[he.edge()] - targetLcr[he.edge()]) /
            targetLcr[he.edge()] *
            (vpg.edgeLengths[kj.edge()] / vpg.edgeLengths[jl.edge()]) *
            (grad_li * vpg.edgeLengths[ik.edge()] -
             grad_ik * vpg.edgeLengths[li.edge()]) /
            vpg.edgeLengths[ik.edge()] / vpg.edgeLengths[ik.edge()];
        // regularizationForce[v] += -P.Kst * localAreaGradient;
      }

      if (P.Ksl != 0) {
        regularizationForce[v] +=
            -P.Ksl * localAreaGradient *
            (vpg.faceAreas[base_he.face()] - targetFaceAreas[base_he.face()]) /
            targetFaceAreas[base_he.face()];
      }

      if (P.Kse != 0) {
        regularizationForce[v] +=
            -P.Kse * edgeGradient *
            (vpg.edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
            targetEdgeLengths[he.edge()];
      }

      /// Patch regularization
      // // the cubic penalty is for regularizing the mesh,
      // // need better physical interpretation or alternative method
      // if (P.Kse != 0) {
      //   double strain =
      //       (vpg.edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
      //       targetEdgeLengths[he.edge()];
      //   regularizationForce[v] +=
      //       -P.Kse * edgeGradient * strain * strain * strain;
      // }
    }
  }
}

} // end namespace ddgsolver