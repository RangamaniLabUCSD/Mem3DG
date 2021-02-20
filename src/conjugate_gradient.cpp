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

void ConjugateGradient::checkParameters() {
  if (f.P.gamma != 0) {
    throw std::runtime_error("gamma has to be 0 for Conjugate Gradient integration!");
  }
}

void ConjugateGradient::status() {
  // map the raw eigen datatype for computation
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);

  // recompute cached values
  f.updateVertexPositions();

  // compute summerized forces
  getForces();

  // compute velocity
  vel_e = f.M * (physicalPressure + DPDPressure) + regularizationForce;

  // compute the L2 error norm
  f.L2ErrorNorm = f.computeL2Norm(vel_e);

  // compute the area contraint error
  dArea = (f.P.Ksg != 0 && !f.mesh->hasBoundary())
              ? abs(f.surfaceArea / f.targetSurfaceArea - 1)
              : 0.0;

  if (f.isReducedVolume) {
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

  // determine conjugate gradient direction, restart after nVertices() cycles
  size_t countCG = 0;
  if (countCG % (f.mesh->nVertices() + 1) == 0) {
    pastNormSq = vel_e.squaredNorm();
    direction = vel_e;
    countCG = 0;
  } else {
    currentNormSq = vel_e.squaredNorm();
    direction = currentNormSq / pastNormSq * direction + vel_e;
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
    backtrack(rho, c1, EXIT, SUCCESS, f.E.potE, vel_e, direction);
  } else {
    pos_e += direction * dt;
    f.time += dt;
  }

  // vertex shift for regularization
  if (f.isVertexShift) {
    f.vertexShift();
  }

  // time stepping on protein density
  if (f.isProtein) {
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
      (f.isReducedVolume ? f.P.Vt : f.P.cam) = VP;
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

void Integrator::backtrack(
    double rho, double c1, bool &EXIT, bool &SUCCESS,
    const double potentialEnergy_pre,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &force,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &direction) {

  // calculate initial energy as reference level
  Eigen::Matrix<double, Eigen::Dynamic, 3> init_position =
      gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  double init_time = f.time;

  // declare variables used in backtracking iterations
  double alpha = dt;
  size_t count = 0;
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);

  pos_e += alpha * direction;
  f.updateVertexPositions();
  f.computeFreeEnergy();

  while (f.E.potE > (potentialEnergy_pre -
                     c1 * alpha * (force.array() * direction.array()).sum())) {
    // while (f.E.potE > potentialEnergy_pre) {
    if (alpha < 1e-12) {
      std::cout << "\nline search failure! Simulation stopped. \n" << std::endl;
      EXIT = true;
      SUCCESS = false;

      // restore entry configuration
      alpha = dt;
      pos_e = init_position;
      f.updateVertexPositions();
      f.computeFreeEnergy();
      f.time = init_time - alpha;

      break;
    }
    alpha *= rho;
    pos_e = init_position + alpha * direction;
    f.updateVertexPositions();
    f.computeFreeEnergy();
    count++;
  }

  if (alpha != dt && verbosity > 3) {
    std::cout << "alpha: " << dt << " -> " << alpha << std::endl;
    std::cout << "L2 norm: " << f.L2ErrorNorm << std::endl;
  }
  f.time = init_time + alpha;
}

void Integrator::pressureConstraintThreshold(bool &EXIT,
                                             const bool isAugmentedLagrangian,
                                             const double dArea,
                                             const double ctol,
                                             double increment) {
  if (f.L2ErrorNorm < tol) {
    if (isAugmentedLagrangian) { // augmented Lagrangian method
      if (dArea < ctol) {        // exit if fulfilled all constraints
        std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[lambdaSG] = [" << f.P.lambdaSG << ", "
                  << "]";
        f.P.lambdaSG += f.P.Ksg * (f.surfaceArea - f.targetSurfaceArea) /
                        f.targetSurfaceArea;
        std::cout << " -> [" << f.P.lambdaSG << "]" << std::endl;
      }
    } else {              // incremental harmonic penalty method
      if (dArea < ctol) { // exit if fulfilled all constraints
        std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[Ksg] = [" << f.P.Ksg << "]";
        f.P.Ksg *= increment;
        std::cout << " -> [" << f.P.Ksg << "]" << std::endl;
      }
    }
  }
}

void Integrator::reducedVolumeThreshold(bool &EXIT,
                                        const bool isAugmentedLagrangian,
                                        const double dArea,
                                        const double dVolume, const double ctol,
                                        double increment) {
  if (f.L2ErrorNorm < tol) {
    if (isAugmentedLagrangian) {            // augmented Lagrangian method
      if (dArea < ctol && dVolume < ctol) { // exit if fulfilled all constraints
        std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[lambdaSG, lambdaV] = [" << f.P.lambdaSG << ", "
                  << f.P.lambdaV << "]";
        f.P.lambdaSG += f.P.Ksg * (f.surfaceArea - f.targetSurfaceArea) /
                        f.targetSurfaceArea;
        f.P.lambdaV +=
            f.P.Kv * (f.volume - f.refVolume * f.P.Vt) / (f.refVolume * f.P.Vt);
        std::cout << " -> [" << f.P.lambdaSG << ", " << f.P.lambdaV << "]"
                  << std::endl;
      }
    } else { // incremental harmonic penalty method
      if (dArea < ctol && dVolume < ctol) { // exit if fulfilled all constraints
        std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
        EXIT = true;
      }

      // iterate if not
      if (dArea > ctol) {
        std::cout << "\n[Ksg] = [" << f.P.Ksg << "]";
        f.P.Ksg *= 1.3;
        std::cout << " -> [" << f.P.Ksg << "]" << std::endl;
      }
      if (dVolume > ctol) {
        std::cout << "\n[Kv] = [" << f.P.Kv << "]";
        f.P.Kv *= 1.3;
        std::cout << " -> [" << f.P.Kv << "]" << std::endl;
      }
    }
  }
}

} // namespace mem3dg
