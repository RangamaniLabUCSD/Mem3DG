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

#include "mem3dg/solver/force.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/meshops.h"

#ifdef MEM3DG_WITH_NETCDF
#include "mem3dg/solver/trajfile.h"
#endif

namespace ddgsolver {
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
    const std::tuple<double, double, double, double, double, double, double,
                     double>
        energy,
    const size_t &verbosity);
#endif

void velocityVerlet(System &f, double dt, double total_time, double tolerance,
                    double closeZone, double increment, double maxKv,
                    double maxKsg, double tSave, double tMollify,
                    const size_t verbosity, std::string inputMesh,
                    std::string outputDir, double init_time,
                    double errorJumpLim) {

  // print out a txt file listing all parameters used
  if (verbosity > 2) {
    getParameterLog(f, dt, total_time, tolerance, tSave, inputMesh, outputDir);
  }

  // initialize variables used in time integration
  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure, newTotalPressure,
      regularizationForce, physicalPressure, DPDForce;

  const double hdt = 0.5 * dt, hdt2 = hdt * dt;
  double totalEnergy, BE, sE, pE, kE, cE, lE, exE, L2ErrorNorm, dArea, dVolume,
      time = init_time; // double dRef;

  size_t nMollify = size_t(tMollify / tSave), frame = 0;

  bool EXIT = false;

  // map the raw eigen datatype for computation
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

// initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  TrajFile fd;
  if (verbosity > 0) {
    fd.createNewFile(outputDir + "/traj.nc", f.mesh, f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(f.mask.cast<int>());
  }
#endif

  // time integration loop
  while (time <= total_time) {

    // compute summerized forces
    getForces(f, physicalPressure, DPDForce, regularizationForce);
    totalPressure.resize(f.mesh.nVertices(), 3);
    totalPressure.setZero();
    newTotalPressure = physicalPressure + DPDForce;

    // measure the error norm and constraint, exit if smaller than tolerance or
    // reach time limit
    L2ErrorNorm = getL2ErrorNorm(physicalPressure);
    dArea = (f.P.Ksg != 0 && !f.mesh.hasBoundary())
                ? abs(f.surfaceArea / f.targetSurfaceArea - 1)
                : 0.0;
    dVolume = (f.P.Kv != 0 && !f.mesh.hasBoundary())
                  ? abs(f.volume / f.refVolume / f.P.Vt - 1)
                  : 0.0;
    if (L2ErrorNorm < tolerance) {
      EXIT = true;
    }

    // compute the free energy of the system
    std::tie(totalEnergy, BE, sE, pE, kE, cE, lE, exE) = getFreeEnergy(f);

    // Save files every tSave period and print some info
    static double lastSave;
    if (time - lastSave >= tSave - 1e-12 || time == total_time ||
        time == init_time || EXIT) {
      lastSave = time;

      // save variable to richData
      if (verbosity > 2) {
        saveRichData(f, physicalPressure, verbosity);
      }

#ifdef MEM3DG_WITH_NETCDF
      // save variable to netcdf traj file
      if (verbosity > 0) {
        saveNetcdfData(f, frame, time, fd, physicalPressure,
                       std::tie(totalEnergy, BE, sE, pE, kE, cE, lE, exE),
                       verbosity);
      }
#endif

      if (verbosity > 1) {
        std::cout << "\n"
                  << "Time: " << time << "\n"
                  << "Frame: " << frame << "\n"
                  << "dArea: " << dArea << "\n"
                  << "dVolume:  " << dVolume << "\n"
                  << "Total energy (exclude V^ext): " << totalEnergy << "\n"
                  << "L2 error norm: " << L2ErrorNorm << "\n"
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

      if (verbosity > 2) {
        char buffer[50];
        sprintf(buffer, "/t=%d", int(time * 100));
        f.richData.write(outputDir + buffer + ".ply");
      }
    }
    // break loop if EXIT flag is on
    if (EXIT) {
      if (verbosity > 0) {
        std::cout << "\n"
                  << "Simulation finished, and data saved to " + outputDir
                  << std::endl;
      }
      break;
    }

    // time stepping on vertex position
    pos_e += vel_e * dt + hdt2 * totalPressure;
    pos_e += regularizationForce * dt;
    vel_e += (totalPressure + newTotalPressure) * hdt;
    totalPressure = newTotalPressure;
    time += dt;
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
}
} // namespace integration
} // namespace ddgsolver
