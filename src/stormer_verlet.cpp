// Copyright (c) 2020 ctlee
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/system.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <iostream>

namespace mem3dg {
namespace integration {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void stormerVerlet(System &f, double dt, double total_time, double tolerance) {
  gcs::FaceData<size_t> faceInd = f.vpg->faceIndices;
  gc::Vector3 totalForce;
  for (size_t i = 0; i < total_time / dt; i++) {
    /*polyscope::registerSurfaceMesh("myMesh",
            ptrVpg->inputVertexPositions,
            ptrMesh->getFaceVertexList());*/
    // polyscope::show();
    gc::EigenMap<double, 3>(f.vel) =
        (gc::EigenMap<double, 3>(f.vpg->inputVertexPositions) -
         gc::EigenMap<double, 3>(f.pastPositions)) /
        dt;
    f.computeBendingPressure();
    f.computeCapillaryPressure();
    f.computeInsidePressure();
    f.computeDPDForces();
    f.computeExternalPressure();

    gcs::VertexData<gc::Vector3> temp = f.vpg->inputVertexPositions;
    for (gcs::Vertex v : f.mesh->vertices()) {
      bool flag = true;
      for (gcs::Face f : v.adjacentFaces()) {
        if (faceInd[f] == 0) {
          flag = true; // change it to false to have time integration excluded
                       // for a facet of vertices
        }
      }
      if (flag == true) {
        f.vpg->inputVertexPositions[v] *= 2;
        totalForce = f.bendingPressure[v] + f.capillaryPressure[v] +
                     f.insidePressure * f.vpg->vertexNormals[v] +
                     f.externalPressure[v] +
                     ((f.dampingForce[v] + f.stochasticForce[v] +
                       f.regularizationForce[v]) /
                      f.vpg->vertexDualAreas[v]);
        f.vpg->inputVertexPositions[v] +=
            totalForce * dt * dt - f.pastPositions[v];
      }
    }
    // std::cout << "total force:  " << totalForce.norm() << std::endl;
    f.updateVertexPositions();
    f.pastPositions = temp;
    f.computeFreeEnergy();

    // std::cout << "energy: " << totalEnergy << std::endl;
    // std::cout << "process: " << int(double(i) / (total_time / dt) * 100) <<
    // "%" << std::endl;
  }
}
} // namespace integration
} // namespace mem3dg
