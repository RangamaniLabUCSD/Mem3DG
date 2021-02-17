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

#ifdef __linux__
#include <sys/time.h>
#endif

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
namespace gc = ::geometrycentral;

void Integrator::velocityVerlet() {
  signal(SIGINT, signalHandler);

  if (verbosity > 1) {
    std::cout << "Running Velocity Verlet integrator ..." << std::endl;
  }

  // check the validity of parameter
  checkParameters("velocity verlet");

  // initialize variables used in time integration
  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure, newTotalPressure;
  const double hdt = 0.5 * dt, hdt2 = hdt * dt;
  double dArea, dVP, init_time = f.time; // double dRef;
  size_t frame = 0;
  bool EXIT = false;

  // initialize variables used if adopting adaptive time step based on mesh size
  double dt_size2_ratio;
  if (isAdaptiveStep) {
    dt_size2_ratio = dt / f.vpg->edgeLengths.raw().minCoeff() /
                     f.vpg->edgeLengths.raw().minCoeff();
  }

// start the timer
#ifdef __linux__
  struct timeval start;
  gettimeofday(&start, NULL);
#endif

  // map the raw eigen datatype for computation
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);

// initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  TrajFile fd;
  if (verbosity > 0) {
    fd.createNewFile(outputDir + trajFileName, *f.mesh, *f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(f.mask.raw().cast<int>());
    fd.writeRefVolume(f.refVolume);
    fd.writeRefSurfArea(f.targetSurfaceArea);
  }
#endif

  // time integration loop
  for (;;) {
    // compute summerized forces
    getForces();
    totalPressure.resize(f.mesh->nVertices(), 3);
    totalPressure.setZero();
    newTotalPressure = physicalPressure + DPDPressure;

    // compute the L2 error norm
    f.L2ErrorNorm =
        f.computeL2Norm(f.M * physicalPressure + regularizationForce);

    // compute the area contraint error
    dArea = (f.P.Ksg != 0 && !f.mesh->hasBoundary())
                ? abs(f.surfaceArea / f.targetSurfaceArea - 1)
                : 0.0;

    if (f.isReducedVolume) {
      // compute volume constraint error
      dVP = (f.P.Kv != 0 && !f.mesh->hasBoundary())
                ? abs(f.volume / f.refVolume / f.P.Vt - 1)
                : 0.0;
    } else {
      // compute pressure constraint error
      dVP = (!f.mesh->hasBoundary()) ? abs(1.0 / f.volume / f.P.cam - 1) : 1.0;
    }

    // exit if under error tol
    if (f.L2ErrorNorm < tol) {
      std::cout << "\nL2 error norm smaller than tol." << std::endl;
      EXIT = true;
    }

    // exit if reached time
    if (f.time > total_time) {
      std::cout << "\nReached time." << std::endl;
      EXIT = true;
    }

    // compute the free energy of the system
    f.computeFreeEnergy();

    // Save files every tSave period and print some info
    static double lastSave;
    if (f.time - lastSave >= tSave - 1e-12 || f.time == init_time || EXIT) {
      lastSave = f.time;

      // save variable to richData and save ply file
      if (verbosity > 3) {
        saveRichData();
        char buffer[50];
        sprintf(buffer, "/frame%d", (int)frame);
        f.richData.write(outputDir + buffer + ".ply");
      }

      // save variable to netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
      if (verbosity > 0) {
        saveNetcdfData(frame, fd);
      }
#endif

      // print in-progress information in the console
      if (verbosity > 1) {
        std::cout << "\n"
                  << "t: " << f.time << ", "
                  << "n: " << frame << "\n"
                  << "dA: " << dArea << ", "
                  << "dVP: " << dVP << ", "
                  << "h: " << abs(f.vpg->inputVertexPositions[f.theVertex].z)
                  << "\n"
                  << "E_total: " << f.E.totalE << "\n"
                  << "|e|L2: " << f.L2ErrorNorm << "\n"
                  << "H: [" << f.H.raw().minCoeff() << ","
                  << f.H.raw().maxCoeff() << "]"
                  << "\n"
                  << "K: [" << f.K.raw().minCoeff() << ","
                  << f.K.raw().maxCoeff() << "]" << std::endl;
        // << "COM: "
        // << gc::EigenMap<double,
        // 3>(f.vpg->inputVertexPositions).colwise().sum() /
        //         f.vpg->inputVertexPositions.raw().rows()
        // << "\n"
      }
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      if (verbosity > 0) {
        std::cout << "Simulation finished, and data saved to " + outputDir
                  << std::endl;
        if (verbosity > 2) {
          saveRichData();
          f.richData.write(outputDir + "/out.ply");
        }
      }
      break;
    }

    // adjust time step if adopt adaptive time step based on mesh size
    if (isAdaptiveStep) {
      double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
      dt = dt_size2_ratio * minMeshLength * minMeshLength;
    }

    // time stepping on vertex position
    pos_e += vel_e * dt + hdt2 * totalPressure;
    pos_e += regularizationForce * dt;
    vel_e += (totalPressure + newTotalPressure) * hdt;
    totalPressure = newTotalPressure;
    f.time += dt;

    // vertex shift for regularization
    if (f.isVertexShift) {
      f.vertexShift();
    }

    // time stepping on protein density
    if (f.isProtein) {
      f.proteinDensity.raw() += -f.P.Bc * f.chemicalPotential.raw() * dt;
    }

    // recompute cached values
    f.updateVertexPositions();

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

void Integrator::getForces() {
  f.computeAllForces();

  physicalPressure = rowwiseScaling(
      f.mask.raw().cast<double>(),
      gc::EigenMap<double, 3>(f.bendingPressure) +
          gc::EigenMap<double, 3>(f.capillaryPressure) +
          f.insidePressure * gc::EigenMap<double, 3>(f.vpg->vertexNormals) +
          gc::EigenMap<double, 3>(f.externalPressure) +
          gc::EigenMap<double, 3>(f.lineTensionPressure));

  regularizationForce =
      rowwiseScaling(f.mask.raw().cast<double>(),
                     gc::EigenMap<double, 3>(f.regularizationForce));

  DPDPressure =
      rowwiseScaling(f.mask.raw().cast<double>(),
                     f.M_inv * (EigenMap<double, 3>(f.dampingForce) +
                                gc::EigenMap<double, 3>(f.stochasticForce)));

  if (!f.mesh->hasBoundary()) {
    removeTranslation(physicalPressure);
    removeRotation(EigenMap<double, 3>(f.vpg->inputVertexPositions),
                   physicalPressure);
    // removeTranslation(DPDPressure);
    // removeRotation(EigenMap<double, 3>(f.vpg->inputVertexPositions),
    // DPDPressure);
  }
}

} // namespace mem3dg
