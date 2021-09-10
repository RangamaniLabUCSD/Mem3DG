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
    std::cout << "Initialized NetCDF file at " << outputDirectory + "/" + trajFileName
              << std::endl;
  }
#endif

  gcs::FaceData<std::size_t> faceInd = system.vpg->faceIndices;
  gc::Vector3 totalForce;
  for (std::size_t i = 0; i < totalTime / characteristicTimeStep; i++) {
    /*polyscope::registerSurfaceMesh("myMesh",
            ptrVpg->inputVertexPositions,
            ptrMesh->getFaceVertexList());*/
    // polyscope::show();
    timeStep = characteristicTimeStep;
    gc::EigenMap<double, 3>(system.velocity) =
        (gc::EigenMap<double, 3>(system.vpg->inputVertexPositions) -
         gc::EigenMap<double, 3>(pastPositions)) /
        timeStep;
    system.computeBendingForce();
    system.computeCapillaryForce();
    system.computeOsmoticForce();
    system.computeDPDForces(timeStep);
    system.computeExternalForce();

    gcs::VertexData<gc::Vector3> temp = system.vpg->inputVertexPositions;
    for (gcs::Vertex v : system.mesh->vertices()) {
      bool flag = true;
      for (gcs::Face f : v.adjacentFaces()) {
        if (faceInd[f] == 0) {
          flag = true; // change it to false to have time integration excluded
                       // for a facet of vertices
        }
      }
      if (flag == true) {
        system.vpg->inputVertexPositions[v] *= 2;
        totalForce = system.forces.mechanicalForceVec[v] + system.forces.dampingForce[v] +
                     system.forces.stochasticForce[v] +
                     system.forces.regularizationForce[v];
        system.vpg->inputVertexPositions[v] +=
            totalForce * timeStep * timeStep - pastPositions[v];
      }
    }
    // std::cout << "total force:  " << totalForce.norm() << std::endl;
    // process the mesh with regularization or mutation
    system.mutateMesh();
    system.updateVertexPositions();
    pastPositions = temp;
    // initialize/update the vertex position of the last
    // iteration
    pastPositions = system.vpg->inputVertexPositions;
    system.computeTotalEnergy();

    // std::cout << "energy: " << totalEnergy << std::endl;
    // std::cout << "process: " << int(double(i) / (total_time / dt) * 100) <<
    // "%" << std::endl;
  }

  return true;
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
