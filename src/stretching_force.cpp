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

template <typename T>
void log(gcs::FaceData<T> face_a, gcs::HalfedgeMesh &mesh, std::string name) {
  for (gcs::Face f : mesh.faces()) {
    std::cout << name << face_a[f] << std::endl;
  }
}

void Force::getStretchingForces() {

  capillaryPressure.fill({0.0, 0.0, 0.0});

  // Alias & gc::EigenMap
  const gcs::FaceData<gc::Vector3> &face_n = vpg.faceNormals;
  const gcs::FaceData<double> &face_a = vpg.faceAreas;
  const gcs::VertexData<size_t> &v_ind = vpg.vertexIndices;
  auto &faceArea_e = vpg.faceAreas.raw();

  surfaceArea = faceArea_e.sum();

  for (gcs::Vertex v : mesh.vertices()) {
    gc::Vector3 localForce{0.0, 0.0, 0.0};
    gc::Vector3 globalForce{0.0, 0.0, 0.0};
    gc::Vector3 edgeForce{0.0, 0.0, 0.0};

    for (gcs::Halfedge he : v.outgoingHalfedges()) {

      gc::Vector3 edgeGradient = -vecFromHalfedge(he, vpg).normalize();
      gcs::Halfedge base_he = he.next();
      gc::Vector3 base_vec = vecFromHalfedge(base_he, vpg);
      gc::Vector3 localAreaGradient = -gc::cross(base_vec, face_n[he.face()]);
      assert((gc::dot(localAreaGradient, vecFromHalfedge(he, vpg))) < 0);

      if (P.Kst != 0) {
        regularizationForce[v] += -P.Kst * localAreaGradient;
      }

      if (P.Ksl != 0) {
        localForce +=
            - P.Ksl * localAreaGradient *
            (face_a[base_he.face()] - targetFaceAreas[base_he.face()]) /
            targetFaceAreas[base_he.face()];
      }

      if (P.Ksg != 0) {
        globalForce += - P.Ksg * localAreaGradient *
                       (surfaceArea - targetSurfaceArea) / targetSurfaceArea;
      }

      if (P.Kse != 0) {
        edgeForce +=
            - P.Kse * edgeGradient *
            (vpg.edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
            targetEdgeLengths[he.edge()];
      }
    }
    capillaryPressure[v] = globalForce / vpg.vertexDualAreas[v];
    regularizationForce[v] = localForce + edgeForce;
  }
}
} // end namespace ddgsolver
