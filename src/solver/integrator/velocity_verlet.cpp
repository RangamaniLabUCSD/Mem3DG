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

#include <Eigen/Core>
#include <csignal>
#include <iostream>
#include <math.h>
#include <pcg_random.hpp>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/integrator/velocity_verlet.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;

bool VelocityVerlet::integrate() {
  signal(SIGINT, signalHandler);

#ifdef __linux__
  // start the timer
  struct timeval start;
  gettimeofday(&start, NULL);
#endif

  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  createNetcdfFile();
  // print to console
  std::cout << "Initialized NetCDF file at " << outputDir + "/" + trajFileName
            << std::endl;
#endif

  // time integration loop
  for (;;) {

    // Evaluate and threhold status data
    status();

    // Save files every tSave period and print some info
    static double lastSave;
    if (f.time - lastSave >= tSave || f.time == init_time || EXIT) {
      lastSave = f.time;
      saveData();
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // step forward
    march();
  }

  // return if physical simulation is sucessful
  if (!SUCCESS) {
    markFileName("_failed");
  }

  // stop the timer and report time spent
#ifdef __linux__
  double duration = getDuration(start);
  if (verbosity > 0) {
    std::cout << "\nTotal integration time: " << duration << " seconds"
              << std::endl;
  }
#endif

  return SUCCESS;
}

void VelocityVerlet::checkParameters() {
  f.meshProcessor.meshMutator.summarizeStatus();
  if (f.meshProcessor.isMeshMutate) {
    mem3dg_runtime_error(
        "Mesh mutations are currently not supported for Velocity Verlet!");
  }
}

void VelocityVerlet::status() {
  // recompute cached values
  f.updateVertexPositions();

  // alias vpg quantities, which should follow the update
  auto physicalForceVec = toMatrix(f.forces.mechanicalForceVec);
  auto physicalForce = toMatrix(f.forces.mechanicalForce);
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(f.vpg->vertexNormals);

  // compute summerized forces
  getForces();

  // Compute total pressure
  // newTotalPressure = rowwiseScalarProduct((physicalForce + DPDForce).array()
  // /
  //                                       f.vpg->vertexDualAreas.raw().array(),
  //                                   vertexAngleNormal_e);
  newTotalPressure =
      rowwiseScalarProduct(DPDForce.array() /
                               f.vpg->vertexDualAreas.raw().array(),
                           vertexAngleNormal_e) +
      (physicalForceVec.array().colwise() /
       f.vpg->vertexDualAreas.raw().array())
          .matrix();

  // compute the area contraint error
  dArea = (f.parameters.tension.Ksg != 0)
              ? abs(f.surfaceArea / f.refSurfaceArea - 1)
              : 0.0;

  if (f.parameters.osmotic.isPreferredVolume) {
    // compute volume constraint error
    dVP = (f.parameters.osmotic.Kv != 0)
              ? abs(f.volume / f.parameters.osmotic.Vt - 1)
              : 0.0;
  } else {
    // compute pressure constraint error
    dVP =
        (!f.mesh->hasBoundary())
            ? abs(f.parameters.osmotic.n / f.volume / f.parameters.osmotic.cam -
                  1.0)
            : 1.0;
  }

  // exit if under error tol
  if (f.mechErrorNorm < tol && f.chemErrorNorm < tol) {
    std::cout << "\nError norm smaller than tol." << std::endl;
    EXIT = true;
  }

  // exit if reached time
  if (f.time > total_time) {
    std::cout << "\nReached time." << std::endl;
    EXIT = true;
  }

  // compute the free energy of the system
  f.computeFreeEnergy();

  // backtracking for error
  finitenessErrorBacktrack();
  if (f.energy.totalE > 1.05 * totalEnergy) {
    std::cout
        << "\nVelocity Verlet: increasing system energy, simulation stopped!"
        << std::endl;
    SUCCESS = false;
  }
}

void VelocityVerlet::march() {

  // map the raw eigen datatype for computation
  auto vel_e = gc::EigenMap<double, 3>(f.velocity);
  auto vel_protein_e = toMatrix(f.proteinVelocity);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);

  // adjust time step if adopt adaptive time step based on mesh size
  if (isAdaptiveStep) {
    double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
    dt = dt_size2_ratio * minMeshLength * minMeshLength;
  }
  double hdt = 0.5 * dt, hdt2 = hdt * dt;
  // f.P.sigma =
  //     sqrt(2 * f.P.gamma * mem3dg::constants::kBoltzmann * f.P.temp / dt);

  // time stepping on vertex position
  previousE = f.energy;
  pos_e += vel_e * dt + hdt2 * totalPressure;
  vel_e += (totalPressure + newTotalPressure) * hdt;
  totalPressure = newTotalPressure;
  f.time += dt;

  // time stepping on protein density
  if (f.parameters.variation.isProteinVariation) {
    vel_protein_e = f.parameters.Bc * f.forces.chemicalPotential.raw();
    f.proteinDensity.raw() += vel_protein_e * dt;
  }

  // process the mesh with regularization or mutation
  f.processMesh();
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
