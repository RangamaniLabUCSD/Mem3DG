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

#ifdef MEM3DG_WITH_NETCDF
#include "mem3dg/solver/trajfile.h"
#endif

namespace mem3dg {
namespace integration {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void getForces(System &f,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &DPDForce,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &regularizationForce);

void saveRichData(
    const System &f,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
    const size_t verbosity);

#ifdef MEM3DG_WITH_NETCDF
void saveNetcdfData(
    const System &f, size_t &frame, const double &time, TrajFile &fd,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
    const size_t &verbosity);
#endif

/**
 * @brief Backtracking algorithm that dynamically adjust step size based on
 * energy evaluation
 * @param f, force object
 * @param dt, initial step size
 * @param rho, discount factor
 * @param c1, constant for Wolfe condtion, between 0 to 1, usually ~ 1e-4
 * @param time, simulation time
 * @param EXIT, exit flag for integration loop
 * @param verbosity, verbosity setting
 * @param potentialEnergy_pre, previous energy evaluation
 * @param force, gradient of the energy
 * @param direction, direction, most likely some function of gradient
 * @return
 */
void backtrack(System &f, const double dt, double rho, double c1, double &time,
               bool &EXIT, const size_t verbosity,
               const double potentialEnergy_pre,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &force,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &direction) {

  // calculate initial energy as reference level
  Eigen::Matrix<double, Eigen::Dynamic, 3> init_position =
      gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);
  double init_time = time;

  // declare variables used in backtracking iterations
  double alpha = dt;
  size_t count = 0;
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  pos_e += alpha * direction;
  f.update_Vertex_positions();
  f.getFreeEnergy();

  while (f.E.potE > (potentialEnergy_pre -
                     c1 * alpha * (force.array() * direction.array()).sum())) {
    // while (f.E.potE > potentialEnergy_pre) {
    if (alpha < 1e-12) {
      std::cout << "\nline search failure! Simulation stopped. \n" << std::endl;
      EXIT = true;

      // restore entry configuration
      alpha = dt;
      pos_e = init_position;
      f.update_Vertex_positions();
      f.getFreeEnergy();
      time = init_time - alpha;

      break;
    }
    alpha *= rho;
    pos_e = init_position + alpha * direction;
    f.update_Vertex_positions();
    f.getFreeEnergy();
    count++;
  }

  if (alpha != dt && verbosity > 3) {
    std::cout << "alpha: " << dt << " -> " << alpha << std::endl;
    std::cout << "L2 norm: " << f.L2ErrorNorm << std::endl;
  }
  time = init_time + alpha;
}

void conjugateGradient(System &f, double dt, double init_time,
                       double total_time, double tSave, double tol, double ctol,
                       const size_t verbosity, std::string outputDir,
                       const bool isBacktrack, const double rho,
                       const double c1, const bool isAugmentedLagrangian,
                       const std::string trajFileName) {

  // initialize variables used in time integration
  Eigen::Matrix<double, Eigen::Dynamic, 3> regularizationForce,
      physicalPressure, DPDForce, direction;
  double dArea, dVolume, currentNormSq, pastNormSq, time = init_time;
  size_t frame = 0;
  bool EXIT = false;

// start the timer
#ifdef __linux__
  struct timeval start;
  gettimeofday(&start, NULL);
#endif
  // map the raw eigen datatype for computation
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  TrajFile fd;
  if (verbosity > 0) {
    fd.createNewFile(outputDir + trajFileName, f.mesh, f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(f.mask.cast<int>());
  }
#endif

  // time integration loop
  for (;;) {

    // compute summerized forces
    getForces(f, physicalPressure, DPDForce, regularizationForce);
    vel_e = physicalPressure + DPDForce + regularizationForce;

    // measure the error norm and constraints, exit if smaller than tolerance
    dArea = (f.P.Ksg != 0 && !f.mesh.hasBoundary())
                ? abs(f.surfaceArea / f.targetSurfaceArea - 1)
                : 0.0;
    dVolume = (f.P.Kv != 0 && !f.mesh.hasBoundary())
                  ? abs(f.volume / f.refVolume / f.P.Vt - 1)
                  : 0.0;
    f.getL2ErrorNorm(physicalPressure);
    if (f.L2ErrorNorm < tol) {
      if (isAugmentedLagrangian) {
        // augmented Lagrangian method
        if (dArea < ctol && dVolume < ctol) {
          std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
          EXIT = true;
        } else {
          std::cout << "\n[lambdaSG, lambdaV] = [" << f.P.lambdaSG << ", "
                    << f.P.lambdaV << "]";
          f.P.lambdaSG += f.P.Ksg * (f.surfaceArea - f.targetSurfaceArea) /
                          f.targetSurfaceArea;
          f.P.lambdaV += f.P.Kv * (f.volume - f.refVolume * f.P.Vt) /
                         (f.refVolume * f.P.Vt);
          std::cout << " -> [" << f.P.lambdaSG << ", " << f.P.lambdaV << "]"
                    << std::endl;
        }
      } else {
        // incremental harmonic penalty method
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
        if (dArea < ctol && dVolume < ctol) {
          std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
          EXIT = true;
        }
      }
    }
    if (time > total_time) {
      std::cout << "\nReached time." << std::endl;
      EXIT = true;
    }

    // compute the free energy of the system
    f.getFreeEnergy();

    // Save files every tSave period and print some info
    static double lastSave;
    if (time - lastSave >= tSave - 1e-12 || time == init_time || EXIT) {
      lastSave = time;

      // save variable to richData and save ply file
      if (verbosity > 3) {
        saveRichData(f, physicalPressure, verbosity);
        char buffer[50];
        sprintf(buffer, "/frame%d", (int)frame);
        f.richData.write(outputDir + buffer + ".ply");
      }

#ifdef MEM3DG_WITH_NETCDF
      // save variable to netcdf traj file
      if (verbosity > 0) {
        saveNetcdfData(f, frame, time, fd, physicalPressure, verbosity);
      }
#endif

      // print in-progress information in the console
      if (verbosity > 1) {
        std::cout << "\n"
                  << "Time: " << time << "\n"
                  << "Frame: " << frame << "\n"
                  << "dArea: " << dArea << "\n"
                  << "dVolume:  " << dVolume << "\n"
                  << "Potential energy (exclude V^ext): " << f.E.potE << "\n"
                  << "L2 error norm: " << f.L2ErrorNorm << "\n"
                  << "COM: "
                  << gc::EigenMap<double, 3>(f.vpg.inputVertexPositions)
                             .colwise()
                             .sum() /
                         f.vpg.inputVertexPositions.raw().rows()
                  << "\n"
                  << "Height: "
                  << abs(f.vpg.inputVertexPositions[f.mesh.vertex(f.ptInd)].z)
                  << "\n";
      }
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      if (verbosity > 0) {
        std::cout << "Simulation finished, and data saved to " + outputDir
                  << std::endl;
        if (verbosity > 2) {
          saveRichData(f, physicalPressure, verbosity);
          f.richData.write(outputDir + "/out.ply");
        }
      }
      break;
    }

    // time stepping on vertex position
    size_t countCG = 0;
    if (countCG % (f.mesh.nVertices() + 1) == 0) {
      pastNormSq = vel_e.squaredNorm();
      direction = vel_e;
      countCG = 0;
    } else {
      currentNormSq = vel_e.squaredNorm();
      direction = currentNormSq / pastNormSq * direction + vel_e;
      pastNormSq = currentNormSq;
      countCG++;
    }
    if (isBacktrack) {
      backtrack(f, dt, rho, c1, time, EXIT, verbosity, f.E.potE, vel_e,
                direction);
    } else {
      pos_e += direction * dt;
      time += dt;
    }
    if (f.isVertexShift) {
      vertexShift(f.mesh, f.vpg, f.mask);
    }

    // time stepping on protein density
    if (f.isProtein) {
      f.proteinDensity.raw() += -f.P.Bc * f.chemicalPotential.raw() * dt;
    }

    // recompute cached values
    f.update_Vertex_positions();

  } // integration

  // stop the timer and report time spent
#ifdef __linux__
  double duration = getDuration(start);
  if (verbosity > 0) {
    std::cout << "\nTotal integration time: " << duration << " seconds"
              << std::endl;
  }
#endif
}

void feedForwardSweep(System &f, std::vector<double> H_, std::vector<double> V_,
                      double dt, double maxTime, double tSave, double tol,
                      double ctol, std::string outputDir,
                      const bool isBacktrack, const double rho, const double c1,
                      const bool isAugmentedLagrangian) {
  const double KV = f.P.Kv, KSG = f.P.Ksg;
  double init_time = 0;
  size_t verbosity = 2;

#ifdef __linux__
  // start the timer
  struct timeval start;
  gettimeofday(&start, NULL);
#endif

  // parameter sweep
  for (double H : H_) {
    for (double V : V_) {
      if (isAugmentedLagrangian) {
        f.P.lambdaSG = 0;
        f.P.lambdaV = 0;
      } else {
        f.P.Kv = KV;
        f.P.Ksg = KSG;
      }
      char buffer[50];
      sprintf(buffer, "/traj_H_%d_V_%d.nc", int(H * 100), int(V * 100));
      f.P.H0 = H;
      f.P.Vt = V;
      std::cout << "\nH0: " << f.P.H0 << std::endl;
      std::cout << "Vt: " << f.P.Vt << std::endl;
      f.update_Vertex_positions();
      conjugateGradient(f, dt, init_time, maxTime, tSave, tol, ctol, verbosity,
                        outputDir, isBacktrack, rho, c1, isAugmentedLagrangian,
                        buffer);
    }
    std::reverse(V_.begin(), V_.end());
  }

  // stop the timer and report time spent
#ifdef __linux__
  double duration = getDuration(start);
  if (verbosity > 0) {
    std::cout << "\nTotal parameter sweep time: " << duration << " seconds"
              << std::endl;
  }
#endif
}

} // namespace integration
} // namespace mem3dg
