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
#include <stdexcept>

#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace gc = ::geometrycentral;

bool BFGS::integrate() {

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
  std::cout << "Initialized integrator and the output trajactory is "
            << outputDir + "/" + trajFileName << std::endl;
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

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // step forward
    march();
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

void BFGS::checkParameters() {
  if (f.P.gamma != 0 || f.P.temp != 0) {
    throw std::runtime_error("DPD has to be turned off for BFGS integration!");
  }
  // if (f.O.isVertexShift) {
  //   throw std::runtime_error(
  //       "Vertex shift is not supported for BFGS integration!");
  // }
  if (!isBacktrack) {
    throw std::runtime_error("Backtracking is required for BFGS integration");
  }
  if (f.P.Bc != 1) {
    throw std::runtime_error(
        "Binding constant should be set to 1 for optimization!");
  }
}

void BFGS::status() {
  // recompute cached values
  f.updateVertexPositions();

  // compute summerized forces
  getForces();

  // update
  if (f.time != init_time || ifRestart) {
    EigenVectorX1D y = -f.F.flatten(physicalForceVec) + pastPhysicalForce;
    double sTy = (s.transpose() * y);
    hess_inv +=
        (s * s.transpose()) * (sTy + y.transpose() * hess_inv * y) / sTy / sTy -
        (hess_inv * y * s.transpose() + s * y.transpose() * hess_inv) / sTy;

    EigenVectorX1D y_protein =
        -f.F.chemicalPotential.raw() + pastPhysicalForce_protein;
    double sTy_protein = (s_protein.transpose() * y_protein);
    hess_inv_protein += (s_protein * s_protein.transpose()) *
                            (sTy_protein + y_protein.transpose() *
                                               hess_inv_protein * y_protein) /
                            sTy_protein / sTy_protein -
                        (hess_inv_protein * y_protein * s_protein.transpose() +
                         s_protein * y_protein.transpose() * hess_inv_protein) /
                            sTy_protein;
  }
  pastPhysicalForce = f.F.flatten(physicalForceVec);
  pastPhysicalForce_protein = f.F.chemicalPotential.raw();
  // std::cout << "if equal: "
  //           << (f.F.unflatten(f.F.flatten(physicalForceVec)).array() ==
  //               physicalForceVec.array())
  //           << std::endl;

  // compute the L1 error norm
  f.L1ErrorNorm = f.computeL1Norm(physicalForce);

  // compute the L1 chemical error norm
  f.L1ChemErrorNorm = f.computeL1Norm(f.F.chemicalPotential.raw());

  // compute the area contraint error
  dArea = (f.P.Ksg != 0) ? abs(f.surfaceArea / f.refSurfaceArea - 1) : 0.0;

  if (f.O.isReducedVolume) {
    // compute volume constraint error
    dVP = (f.P.Kv != 0) ? abs(f.volume / f.refVolume / f.P.Vt - 1) : 0.0;
    // thresholding, exit if fulfilled and iterate if not
    reducedVolumeThreshold(EXIT, isAugmentedLagrangian, dArea, dVP, ctol, 1.3);
  } else {
    // compute pressure constraint error
    dVP = (!f.mesh->hasBoundary()) ? abs(1 / f.volume / f.P.cam - 1) : 1.0;
    // thresholding, exit if fulfilled and iterate if not
    pressureConstraintThreshold(EXIT, isAugmentedLagrangian, dArea, ctol, 1.3);
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

void BFGS::march() {
  if (f.time == lastSave && f.time != init_time) {
    // process the mesh with regularization or mutation
    f.processMesh();
    f.time += 1e-10 * dt;
    ifRestart = true;
    hess_inv.setIdentity();
    hess_inv_protein.setIdentity();
  } else {
    // map the raw eigen datatype for computation
    auto vel_e = f.F.toMatrix(f.vel);
    auto vel_protein_e = f.F.toMatrix(f.vel_protein);

    vel_e = f.F.unflatten(hess_inv * f.F.flatten(physicalForceVec));
    vel_protein_e = hess_inv_protein * f.F.chemicalPotential.raw();

    // adjust time step if adopt adaptive time step based on mesh size
    if (isAdaptiveStep) {
      double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
      dt = dt_size2_ratio * maxForce * minMeshLength * minMeshLength /
           (f.O.isShapeVariation
                ? physicalForce.cwiseAbs().maxCoeff()
                : f.F.chemicalPotential.raw().cwiseAbs().maxCoeff());
    }

    // time stepping on vertex position
    previousE = f.E;
    double alpha = backtrack(f.E.potE, vel_e, vel_protein_e, rho, c1);
    s = alpha * f.F.flatten(vel_e);
    s_protein = alpha * vel_protein_e;

    // regularization
    if ((f.P.Kse != 0) || (f.P.Ksl != 0) || (f.P.Kst != 0)) {
      f.computeRegularizationForce();
      f.vpg->inputVertexPositions.raw() += f.F.regularizationForce.raw();
    }
  }
}

} // namespace mem3dg
