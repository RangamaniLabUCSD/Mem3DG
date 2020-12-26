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
#include <sys/time.h>

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

/**
 * @brief Summerize forces into 3 categories: physcialPressure, DPDForce and
 * regularizationForce. Note that the forces has been removed rigid body mode
 * and masked for integration
 *
 * @param f
 * @param physicalPressure
 * @param DPDForce
 * @param regularizationForce
 * @return
 */
void getForces(System &f,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &DPDForce,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &regularizationForce) {
  f.getAllForces();

  physicalPressure =
      rowwiseScaling(f.mask.cast<double>(),
                     gc::EigenMap<double, 3>(f.bendingPressure) +
                         gc::EigenMap<double, 3>(f.capillaryPressure) +
                         gc::EigenMap<double, 3>(f.insidePressure) +
                         gc::EigenMap<double, 3>(f.externalPressure) +
                         gc::EigenMap<double, 3>(f.lineTensionPressure));

  regularizationForce = rowwiseScaling(
      f.mask.cast<double>(), gc::EigenMap<double, 3>(f.regularizationForce));

  DPDForce =
      rowwiseScaling(f.mask.cast<double>(),
                     f.M_inv * (EigenMap<double, 3>(f.dampingForce) +
                                gc::EigenMap<double, 3>(f.stochasticForce)));

  if (!f.mesh.hasBoundary()) {
    removeTranslation(physicalPressure);
    removeRotation(EigenMap<double, 3>(f.vpg.inputVertexPositions),
                   physicalPressure);
    // removeTranslation(DPDForce);
    // removeRotation(EigenMap<double, 3>(f.vpg.inputVertexPositions),
    // DPDForce);
  }
}

/**
 * @brief Save data to richData
 * @param f, force object
 * @param physcialPressure, physical pressre eigen matrix
 * @param verbosity, verbosity setting
 * @return
 */
void saveRichData(
    const System &f,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
    const size_t verbosity) {
  gcs::VertexData<double> H(f.mesh), H0(f.mesh), fn(f.mesh), f_ext(f.mesh),
      fb(f.mesh), fl(f.mesh), ft(f.mesh);

  H.fromVector(f.H);
  H0.fromVector(f.H0);
  fn.fromVector(rowwiseDotProduct(
      physicalPressure, gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
  f_ext.fromVector(
      rowwiseDotProduct(gc::EigenMap<double, 3>(f.externalPressure),
                        gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
  fb.fromVector(
      rowwiseDotProduct(EigenMap<double, 3>(f.bendingPressure),
                        gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
  fl.fromVector(
      rowwiseDotProduct(EigenMap<double, 3>(f.lineTensionPressure),
                        gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
  ft.fromVector((rowwiseDotProduct(EigenMap<double, 3>(f.capillaryPressure),
                                   gc::EigenMap<double, 3>(f.vpg.vertexNormals))
                     .array() /
                 f.H.array() / 2)
                    .matrix());

  f.richData.addVertexProperty("mean_curvature", H);
  f.richData.addVertexProperty("spon_curvature", H0);
  f.richData.addVertexProperty("external_pressure", f_ext);
  f.richData.addVertexProperty("physical_pressure", fn);
  f.richData.addVertexProperty("capillary_pressure", ft);
  f.richData.addVertexProperty("bending_pressure", fb);
  f.richData.addVertexProperty("line_tension_pressure", fl);
}

#ifdef MEM3DG_WITH_NETCDF
/**
 * @brief Save data to netcdf traj file
 * @param f, force object
 * @param frame, frame index of netcdf traj file
 * @param time, simulation time
 * @param fd, netcdf trajFile object
 * @param physcialPressure, physical pressre eigen matrix
 * @param energy, components of energy - totalE, BE, sE, pE, kE, cE, lE,
 * exE
 * @param verbosity, verbosity setting
 * @return
 */
void saveNetcdfData(
    const System &f, size_t &frame, const double &time, TrajFile &fd,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
    const size_t &verbosity) {

  Eigen::Matrix<double, Eigen::Dynamic, 1> fn, f_ext, fb, fl, ft;

  fn = rowwiseDotProduct(physicalPressure,
                         gc::EigenMap<double, 3>(f.vpg.vertexNormals));
  f_ext = rowwiseDotProduct(gc::EigenMap<double, 3>(f.externalPressure),
                            gc::EigenMap<double, 3>(f.vpg.vertexNormals));
  fb = rowwiseDotProduct(EigenMap<double, 3>(f.bendingPressure),
                         gc::EigenMap<double, 3>(f.vpg.vertexNormals));
  fl = rowwiseDotProduct(EigenMap<double, 3>(f.lineTensionPressure),
                         gc::EigenMap<double, 3>(f.vpg.vertexNormals));
  ft = (rowwiseDotProduct(EigenMap<double, 3>(f.capillaryPressure),
                          gc::EigenMap<double, 3>(f.vpg.vertexNormals))
            .array() /
        f.H.array() / 2)
           .matrix();

  frame = fd.getNextFrameIndex();
  fd.writeTime(frame, time);
  fd.writeCoords(frame, EigenMap<double, 3>(f.vpg.inputVertexPositions));
  fd.writeVelocity(frame, EigenMap<double, 3>(f.vel));
  fd.writeAngles(frame, f.vpg.cornerAngles.raw());
  fd.writeProteinDensity(frame, f.proteinDensity.raw());
  fd.writeMeanCurvature(frame, f.H);
  fd.writeSponCurvature(frame, f.H0);
  fd.writeH_H0_diff(frame,
                    ((f.H - f.H0).array() * (f.H - f.H0).array()).matrix());
  fd.writeExternalPressure(frame, f_ext);
  fd.writePhysicalPressure(frame, fn);
  fd.writeCapillaryPressure(frame, ft);
  fd.writeBendingPressure(frame, fb);
  fd.writeLinePressure(frame, fl);
  fd.writeBendEnergy(frame, f.E.BE);
  fd.writeSurfEnergy(frame, f.E.sE);
  fd.writePressEnergy(frame, f.E.pE);
  fd.writeKineEnergy(frame, f.E.kE);
  fd.writeChemEnergy(frame, f.E.cE);
  fd.writeLineEnergy(frame, f.E.lE);
  fd.writeTotalEnergy(frame, f.E.totalE);
  fd.writeL2ErrorNorm(frame, f.L2ErrorNorm);
}
#endif

void velocityVerlet(System &f, double dt, double init_time, double total_time,
                    double tSave, double tolerance, const size_t verbosity,
                    std::string outputDir) {

  // initialize variables used in time integration
  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure, newTotalPressure,
      regularizationForce, physicalPressure, DPDForce;
  struct timeval start;
  const double hdt = 0.5 * dt, hdt2 = hdt * dt;
  double dArea, dVolume, time = init_time; // double dRef;
  size_t frame = 0;
  bool EXIT = false;

  // start the timer
  gettimeofday(&start, NULL);

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
  for (;;) {
    // compute summerized forces
    getForces(f, physicalPressure, DPDForce, regularizationForce);
    totalPressure.resize(f.mesh.nVertices(), 3);
    totalPressure.setZero();
    newTotalPressure = physicalPressure + DPDForce;

    // measure the error norm and constraint, exit if smaller than tolerance or
    // reach time limit
    f.getL2ErrorNorm(physicalPressure);
    dArea = (f.P.Ksg != 0 && !f.mesh.hasBoundary())
                ? abs(f.surfaceArea / f.targetSurfaceArea - 1)
                : 0.0;
    dVolume = (f.P.Kv != 0 && !f.mesh.hasBoundary())
                  ? abs(f.volume / f.refVolume / f.P.Vt - 1)
                  : 0.0;
    if (f.L2ErrorNorm < tolerance) {
      std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
      EXIT = true;
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

      // save variable to netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
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
                  << "Total energy (exclude V^ext): " << f.E.totalE << "\n"
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

  // stop the timer and report time spent
  double duration = getDuration(start);
  if (verbosity > 0) {
    std::cout << "\nTotal integration time: " << duration << " seconds"
              << std::endl;
  }
}
} // namespace integration
} // namespace ddgsolver
