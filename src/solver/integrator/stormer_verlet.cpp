// Copyright (c) 2020 ctlee
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "mem3dg/solver/integrator/stormer_verlet.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include <iostream>

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

bool StormerVerlet::integrate() {

  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  if (verbosity > 0) {
    // createNetcdfFile();
    createMutableNetcdfFile();
    // print to console
    std::cout << "Initialized NetCDF file at " << outputDir + "/" + trajFileName
              << std::endl;
  }
#endif

  gcs::FaceData<std::size_t> faceInd = f.vpg->faceIndices;
  gc::Vector3 totalForce;
  for (std::size_t i = 0; i < total_time / dt; i++) {
    /*polyscope::registerSurfaceMesh("myMesh",
            ptrVpg->inputVertexPositions,
            ptrMesh->getFaceVertexList());*/
    // polyscope::show();
    gc::EigenMap<double, 3>(f.velocity) =
        (gc::EigenMap<double, 3>(f.vpg->inputVertexPositions) -
         gc::EigenMap<double, 3>(pastPositions)) /
        dt;
    f.computeBendingForce();
    f.computeCapillaryForce();
    f.computeOsmoticForce();
    f.computeDPDForces(dt);
    f.computeExternalForce();

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
        totalForce = f.forces.mechanicalForceVec[v] + f.forces.dampingForce[v] +
                     f.forces.stochasticForce[v] +
                     f.forces.regularizationForce[v];
        f.vpg->inputVertexPositions[v] +=
            totalForce * dt * dt - pastPositions[v];
      }
    }
    // std::cout << "total force:  " << totalForce.norm() << std::endl;
    // process the mesh with regularization or mutation
    f.mutateMesh();
    f.updateVertexPositions();
    pastPositions = temp;
    // initialize/update the vertex position of the last
    // iteration
    pastPositions = f.vpg->inputVertexPositions;
    f.computeFreeEnergy();

    // std::cout << "energy: " << totalEnergy << std::endl;
    // std::cout << "process: " << int(double(i) / (total_time / dt) * 100) <<
    // "%" << std::endl;
  }

  return true;
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
