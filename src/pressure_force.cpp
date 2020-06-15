
#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/IterativeLinearSolvers>

#include "ddgsolver/force.h"
#include "ddgsolver/meshops.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getPressureForces(double &Kv, double &Vt) {
  // pressureForces.fill({0, 0, 0});
  volume = 0;
  double face_volume;
  gcs::FaceData<int> sign_of_volume(mesh);
  for (gcs::Face f : mesh.faces()) {
    face_volume = signedVolumeFromFace(f, vpg);
    volume += face_volume;
    if (face_volume < 0) {
      sign_of_volume[f] = -1;
    } else {
      sign_of_volume[f] = 1;
    }
  }

  std::cout << "total volume:  " << volume / targetVolume / Vt << std::endl;

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
      pressureForces[v] +=
          -0.5 * Kv * (volume - targetVolume * Vt) / (targetVolume * Vt) * dVdx;
      // for (size_t i = 0; i < 3; i++) {
      //	pressureForces(v_ind[v], i) += dVdx[i];
      //}
      // force.row(v_ind[v]) << dVdx.x, dVdx.y, dVdx.z;
    }
  }
}
} // end namespace ddgsolver
