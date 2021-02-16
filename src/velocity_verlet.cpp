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
namespace integration {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

/**
 * @brief Summerize forces into 3 categories: physcialPressure, DPDPressure and
 * regularizationForce. Note that the forces has been removed rigid body mode
 * and masked for integration
 *
 * @param f
 * @param physicalPressure
 * @param DPDPressure
 * @param regularizationForce
 * @return
 */
void getForces(System &f,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &DPDPressure,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &regularizationForce) {
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

/**
 * @brief Save data to richData
 * @param f, force object
 * @param physcialPressure, physical pressre eigen matrix
 * @param verbosity, verbosity setting
 * @return
 */
void saveRichData(
    System &f, const Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
    const size_t verbosity) {
  gcs::VertexData<double> fn(*f.mesh), f_ext(*f.mesh), fb(*f.mesh), fl(*f.mesh),
      ft(*f.mesh);

  fn.fromVector(rowwiseDotProduct(
      physicalPressure, gc::EigenMap<double, 3>(f.vpg->vertexNormals)));
  f_ext.fromVector(
      rowwiseDotProduct(gc::EigenMap<double, 3>(f.externalPressure),
                        gc::EigenMap<double, 3>(f.vpg->vertexNormals)));
  fb.fromVector(
      rowwiseDotProduct(EigenMap<double, 3>(f.bendingPressure),
                        gc::EigenMap<double, 3>(f.vpg->vertexNormals)));
  fl.fromVector(
      rowwiseDotProduct(EigenMap<double, 3>(f.lineTensionPressure),
                        gc::EigenMap<double, 3>(f.vpg->vertexNormals)));
  ft.fromVector(
      (rowwiseDotProduct(EigenMap<double, 3>(f.capillaryPressure),
                         gc::EigenMap<double, 3>(f.vpg->vertexNormals))
           .array() /
       f.H.raw().array() / 2)
          .matrix());

  f.richData.addVertexProperty("mean_curvature", f.H);
  f.richData.addVertexProperty("gauss_curvature", f.K);
  f.richData.addVertexProperty("spon_curvature", f.H0);
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

  Eigen::Matrix<double, Eigen::Dynamic, 1> fn, f_ext, fb, fl, ft, fp;

  f_ext = rowwiseDotProduct(gc::EigenMap<double, 3>(f.externalPressure),
                            gc::EigenMap<double, 3>(f.vpg->vertexNormals));
  fb = rowwiseDotProduct(EigenMap<double, 3>(f.bendingPressure),
                         gc::EigenMap<double, 3>(f.vpg->vertexNormals));
  fl = rowwiseDotProduct(EigenMap<double, 3>(f.lineTensionPressure),
                         gc::EigenMap<double, 3>(f.vpg->vertexNormals));
  ft = (rowwiseDotProduct(EigenMap<double, 3>(f.capillaryPressure),
                          gc::EigenMap<double, 3>(f.vpg->vertexNormals)));
  fp.setConstant(f.mesh->nVertices(), 1, f.insidePressure);
  fn = rowwiseDotProduct(physicalPressure,
                         gc::EigenMap<double, 3>(f.vpg->vertexNormals));

  frame = fd.getNextFrameIndex();

  // write time
  fd.writeTime(frame, time);

  // write geometry
  fd.writeCoords(frame, EigenMap<double, 3>(f.vpg->inputVertexPositions));
  fd.writeVolume(frame, f.volume);
  fd.writeSurfArea(frame, f.surfaceArea);
  fd.writeMeanCurvature(frame, f.H.raw());
  fd.writeGaussCurvature(frame, f.K.raw());
  fd.writeSponCurvature(frame, f.H0.raw());
  fd.writeHeight(frame, abs(f.vpg->inputVertexPositions[f.theVertex].z));
  // fd.writeAngles(frame, f.vpg.cornerAngles.raw());
  // fd.writeH_H0_diff(frame,
  //                   ((f.H - f.H0).array() * (f.H -
  //                   f.H0).array()).matrix());

  // write velocity
  fd.writeVelocity(frame, EigenMap<double, 3>(f.vel));
  if (f.isProtein) {
    fd.writeProteinDensity(frame, f.proteinDensity.raw());
  }

  // write pressures
  fd.writeBendingPressure(frame, fb);
  fd.writeCapillaryPressure(frame, ft);
  fd.writeLinePressure(frame, fl);
  fd.writeInsidePressure(frame, fp);
  fd.writeExternalPressure(frame, f_ext);
  fd.writePhysicalPressure(frame, fn);

  // write energies
  fd.writeBendEnergy(frame, f.E.BE);
  fd.writeSurfEnergy(frame, f.E.sE);
  fd.writePressEnergy(frame, f.E.pE);
  fd.writeKineEnergy(frame, f.E.kE);
  fd.writeChemEnergy(frame, f.E.cE);
  fd.writeLineEnergy(frame, f.E.lE);
  fd.writeTotalEnergy(frame, f.E.totalE);

  // write Norms
  fd.writeL2ErrorNorm(frame, f.L2ErrorNorm);
  fd.writeL2BendNorm(frame,
                     f.computeL2Norm(EigenMap<double, 3>(f.bendingPressure)));
  fd.writeL2SurfNorm(frame,
                     f.computeL2Norm(EigenMap<double, 3>(f.capillaryPressure)));
  fd.writeL2PressNorm(
      frame, f.computeL2Norm(f.insidePressure *
                         gc::EigenMap<double, 3>(f.vpg->vertexNormals)));
}
#endif

void velocityVerlet(System &f, double dt, double total_time, double tSave,
                    double tolerance, const size_t verbosity,
                    const bool isAdaptiveStep, std::string outputDir) {
  signal(SIGINT, signalHandler);

  // initialize variables used in time integration
  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure, newTotalPressure,
      regularizationForce, physicalPressure, DPDPressure;
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
    fd.createNewFile(outputDir + "/traj.nc", *f.mesh, *f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(f.mask.raw().cast<int>());
    fd.writeRefVolume(f.refVolume);
    fd.writeRefSurfArea(f.targetSurfaceArea);
  }
#endif

  // time integration loop
  for (;;) {
    // compute summerized forces
    getForces(f, physicalPressure, DPDPressure, regularizationForce);
    totalPressure.resize(f.mesh->nVertices(), 3);
    totalPressure.setZero();
    newTotalPressure = physicalPressure + DPDPressure;

    // compute the L2 error norm
    f.L2ErrorNorm = f.computeL2Norm(f.M * physicalPressure + regularizationForce);

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
    if (f.L2ErrorNorm < tolerance) {
      std::cout << "\nL2 error norm smaller than tolerance." << std::endl;
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
        saveRichData(f, physicalPressure, verbosity);
        char buffer[50];
        sprintf(buffer, "/frame%d", (int)frame);
        f.richData.write(outputDir + buffer + ".ply");
      }

      // save variable to netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
      if (verbosity > 0) {
        saveNetcdfData(f, frame, f.time, fd, physicalPressure, verbosity);
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
          saveRichData(f, physicalPressure, verbosity);
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
} // namespace integration
} // namespace mem3dg
