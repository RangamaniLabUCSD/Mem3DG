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

#include "Eigen/src/Core/util/Constants.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace gc = ::geometrycentral;

bool ConjugateGradient::integrate() {

  signal(SIGINT, signalHandler);

#ifdef __linux__
  // start the timer
  struct timeval start;
  gettimeofday(&start, NULL);
#endif

  // time integration loop
  for (;;) {

    // Evaluate and threhold status data
    status();

    // Save files every tSave period and print some info
    static double lastSave;
    if (f.time - lastSave >= tSave - 1e-12 || f.time == init_time || EXIT) {
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

void ConjugateGradient::checkParameters() {
  if (f.P.gamma != 0 || f.P.temp != 0) {
    throw std::runtime_error("DPD has to be turned off for CG integration!");
  }
}

void ConjugateGradient::status() {
  // recompute cached values
  f.updateVertexPositions();

  // compute summerized forces
  getForces();

  // compute the L1 error norm
  f.L1ErrorNorm = f.computeL1Norm(physicalForce);

  // compute the area contraint error
  dArea = (f.P.Ksg != 0 && !f.mesh->hasBoundary())
              ? abs(f.surfaceArea / f.targetSurfaceArea - 1)
              : 0.0;

  if (f.O.isReducedVolume) {
    // compute volume constraint error
    dVP = (f.P.Kv != 0 && !f.mesh->hasBoundary())
              ? abs(f.volume / f.refVolume / f.P.Vt - 1)
              : 0.0;
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
}

void ConjugateGradient::march() {

  // map the raw eigen datatype for computation
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(f.vpg->vertexNormals);
  // typedef gc::EigenVectorMap_T<double, 3, Eigen::RowMajor>(*EigenMap3)(gcs::VertexData<gc::Vector3>);
  // EigenMap3 Map3 = gc::EigenMap<double, 3>;

  // determine conjugate gradient direction, restart after nVertices() cycles
  size_t countCG = 0;
  if (countCG % (f.mesh->nVertices() + 1) == 0) {
    pastNormSq = physicalForce.squaredNorm();
    vel_e = rowwiseScaling(physicalForce,
                           gc::EigenMap<double, 3>(f.vpg->vertexNormals));
    countCG = 0;
  } else {
    currentNormSq = physicalForce.squaredNorm();
    vel_e = currentNormSq / pastNormSq * vel_e +
            rowwiseScaling(physicalForce,
                           gc::EigenMap<double, 3>(f.vpg->vertexNormals));
    pastNormSq = currentNormSq;
    countCG++;
  }

  // adjust time step if adopt adaptive time step based on mesh size
  if (isAdaptiveStep) {
    double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
    dt = dt_size2_ratio * minMeshLength * minMeshLength;
  }

  // time stepping on vertex position
  if (isBacktrack) {
    backtrack(rho, c1, EXIT, SUCCESS, f.E.potE, physicalForce, vel_e);
  } else {
    pos_e += vel_e * dt;
    f.time += dt;
  }

  // time stepping on protein density
  if (f.O.isProtein) {
    f.proteinDensity.raw() += -f.P.Bc * f.chemicalPotential.raw() * dt;
  }
}

void FeedForwardSweep::sweep() {
#ifdef __linux__
  // start the timer
  struct timeval start;
  gettimeofday(&start, NULL);
#endif

  // initialize variables
  const double KV = f.P.Kv, KSG = f.P.Ksg, init_time = 0.0;
  const size_t verbosity = 2;

  // initialize variables used if adopting adaptive time step based on mesh size
  double dt_size2_ratio;
  if (isAdaptiveStep) {
    dt_size2_ratio = dt / f.vpg->edgeLengths.raw().minCoeff() /
                     f.vpg->edgeLengths.raw().minCoeff();
  }

  // parameter sweep
  for (double H : H_) {
    for (double VP : VP_) {
      // reset parameters
      if (isAugmentedLagrangian) {
        f.P.lambdaSG = 0;
        f.P.lambdaV = 0;
      } else {
        f.P.Kv = KV;
        f.P.Ksg = KSG;
      }

      // adjust time step if adopt adaptive time step based on mesh size
      if (isAdaptiveStep) {
        double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
        dt = dt_size2_ratio * minMeshLength * minMeshLength;
      }

      // initialize trajectory file name
      char buffer[50];
      sprintf(buffer, "/traj_H_%d_VP_%d.nc", int(H * 100), int(VP * 100));
      trajFileName = buffer;

      // update sweeping paraemters
      f.P.H0 = H;
      (f.O.isReducedVolume ? f.P.Vt : f.P.cam) = VP;
      std::cout << "\nH0: " << f.P.H0 << std::endl;
      std::cout << "VP: " << VP << std::endl;

      // recalculate cached values
      f.updateVertexPositions();

      // rerun CG optimization
      integrate();
    }
    // reverse the order of inner loop to ensure phase space closeness
    std::reverse(VP_.begin(), VP_.end());
  }

#ifdef __linux__
  // stop the timer and report time spent
  double duration = getDuration(start);
  if (verbosity > 0) {
    std::cout << "\nTotal parameter sweep time: " << duration << " seconds"
              << std::endl;
  }
#endif
}

} // namespace mem3dg
