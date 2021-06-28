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
#include <iostream>
#include <math.h>
#include <pcg_random.hpp>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/forward_euler.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;

bool Euler::integrate() {

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
    if (f.time - lastSave >= tSave || f.time == init_time || EXIT) {
      lastSave = f.time;
      saveData();
    }

    // Process mesh every tProcessMesh period
    if (f.time - lastProcessMesh > tProcessMesh) {
      lastProcessMesh = f.time;
      f.processMesh();
      f.updateVertexPositions(false);
    }

    // update geodesics every tUpdateGeodesics period
    if (f.time - lastUpdateGeodesics > tUpdateGeodesics) {
      lastUpdateGeodesics = f.time;
      f.updateVertexPositions(true);
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // step forward
    if (f.time == lastProcessMesh || f.time == lastUpdateGeodesics) {
      f.time += 1e-10 * dt;
    } else {
      march();
    }
  }

  // return if optimization is sucessful
  if (!SUCCESS) {
    if (tol == 0) {
      markFileName("_most");
    } else {
      markFileName("_failed");
    }
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

void Euler::checkParameters() {
  if (f.P.gamma != 0 || f.P.temp != 0) {
    mem3dg_runtime_error("DPD has to be turned off for euler integration!");
  }
  if (isBacktrack) {
    if (rho >= 1 || rho <= 0 || c1 >= 1 || c1 <= 0) {
      mem3dg_runtime_error("To backtrack, 0<rho<1 and 0<c1<1!");
    }
  }
}

void Euler::status() {
  // compute summerized forces
  getForces();

  // compute the area contraint error
  dArea = abs(f.surfaceArea / f.refSurfaceArea - 1);
  dVP = (f.O.isReducedVolume) ? abs(f.volume / f.refVolume / f.P.Vt - 1)
                              : abs(1.0 / f.volume / f.P.cam - 1);

  // exit if under error tolerance
  if (f.mechErrorNorm < tol && f.chemErrorNorm < tol) {
    std::cout << "\nError norm smaller than tolerance." << std::endl;
    EXIT = true;
  }

  // exit if reached time
  if (f.time > total_time) {
    std::cout << "\nReached time." << std::endl;
    EXIT = true;
    SUCCESS = false;
  }

  // compute the free energy of the system
  f.computeFreeEnergy();

  // backtracking for error
  finitenessErrorBacktrack();
}

void Euler::march() {
  // map the raw eigen datatype for computation
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto vel_protein_e = toMatrix(f.vel_protein);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  auto physicalForceVec = toMatrix(f.F.mechanicalForceVec);
  auto physicalForce = toMatrix(f.F.mechanicalForce);

  // compute force, which is equivalent to velocity
  vel_e = physicalForceVec;
  vel_protein_e = f.P.Bc * f.F.chemicalPotential.raw();

  // adjust time step if adopt adaptive time step based on mesh size
  if (isAdaptiveStep) {
    double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
    dt = dt_size2_ratio * maxForce * minMeshLength * minMeshLength /
         (f.O.isShapeVariation ? physicalForce.cwiseAbs().maxCoeff()
                               : vel_protein_e.cwiseAbs().maxCoeff());
  }

  // time stepping on vertex position
  previousE = f.E;
  if (isBacktrack) {
    backtrack(f.E.potE, vel_e, vel_protein_e, rho, c1);
  } else {
    pos_e += vel_e * dt;
    f.proteinDensity.raw() += vel_protein_e * dt;
    f.time += dt;
  }

  // regularization
  if ((f.P.Kse != 0) || (f.P.Ksl != 0) || (f.P.Kst != 0)) {
    f.computeRegularizationForce();
    f.vpg->inputVertexPositions.raw() += f.F.regularizationForce.raw();
  }

  // recompute cached values
  f.updateVertexPositions(false);
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
