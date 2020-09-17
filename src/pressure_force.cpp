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

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/solver/force.h"
#include "mem3dg/solver/meshops.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getPressureForces() {
  insidePressure.fill({0, 0, 0});
  volume = 0;
  double face_volume;
  gcs::FaceData<int> sign_of_volume(mesh);

  for (gcs::Face f : mesh.faces()) {
    face_volume = signedVolumeFromFace(
        f, vpg, refVpg.inputVertexPositions[mesh.vertex(P.ptInd)]);
    volume += face_volume;
    if (face_volume < 0) {
      sign_of_volume[f] = -1;
    } else {
      sign_of_volume[f] = 1;
    }
  }

  for (gcs::Vertex v : mesh.vertices()) {
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gcs::Halfedge base_he = he.next();
      gc::Vector3 p1 = vpg.inputVertexPositions[base_he.vertex()];
      gc::Vector3 p2 = vpg.inputVertexPositions[base_he.next().vertex()];
      gc::Vector3 dVdx = 0.5 * gc::cross(p1, p2) / 3.0;
      assert(gc::dot(dVdx, vpg.inputVertexPositions[v] - p1) *
                 sign_of_volume[he.face()] >
             0);
      // std::cout << "i am here" << (gc::dot(dVdx, vpg.inputVertexPositions[v]
      // - p1) < 0) <<  std::endl;
      // dVdx *= sign_of_volume[he.face()];
      insidePressure[v] +=
          - P.Kv * (volume - refVolume * P.Vt) / (refVolume * P.Vt) * dVdx;
      // for (size_t i = 0; i < 3; i++) {
      //	insidePressure(v_ind[v], i) += dVdx[i];
      //}
      // force.row(v_ind[v]) << dVdx.x, dVdx.y, dVdx.z;
    }
  }
}
} // end namespace ddgsolver
