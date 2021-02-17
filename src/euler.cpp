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

bool Integrator::euler(const bool isBacktrack, const double rho,
                       const double c1) {

  signal(SIGINT, signalHandler);

  if (verbosity > 1) {
    std::cout << "Running Forward Euler (steepest descent) propagator ..."
              << std::endl;
  }

  // check the validity of parameter
  checkParameters("euler");

  // initialize variables used in time integration
  double dArea, dVP, init_time = f.time;
  size_t frame = 0;
  bool EXIT = false, SUCCESS = true;

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
    } else {
      // compute pressure constraint error
      dVP = (!f.mesh->hasBoundary()) ? abs(1.0 / f.volume / f.P.cam - 1) : 1.0;
    }

    // exit if under error tolerance
    if (f.L2ErrorNorm < tol) {
      std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
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

#ifdef MEM3DG_WITH_NETCDF
      // save variable to netcdf traj file
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
        std::cout << "Simulation " << (SUCCESS ? "finished" : "failed")
                  << ", and data saved to " + outputDir << std::endl;
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
    if (isBacktrack) {
      backtrack(rho, c1, EXIT, SUCCESS, f.E.potE, vel_e, vel_e);
    } else {
      pos_e += vel_e * dt;
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

  // return if optimization is sucessful
  if (!SUCCESS) {
    markFileName("_failed");
  }
  return SUCCESS;
}
} // namespace mem3dg
