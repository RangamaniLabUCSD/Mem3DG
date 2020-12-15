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
#include <assert.h>
#include <iostream>
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

/**
 * @brief Summerize forces into 3 categories: physcialPressure, DPDForce and
 * regularizationForce. Note that the forces has been removed rigid body mode
 * and masked for integration
 *
 * @param f
 * @param physicalPressure
 * @param regularizationForce
 * @return
 */
void getForces(Force &f,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &DPDForce,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &regularizationForce) {
  if (f.mesh.hasBoundary()) {
    f.getPatchForces();
  } else {
    f.getVesicleForces();
  }
  f.getDPDForces();
  f.getExternalForces();

  if (f.isProtein) {
    f.getChemicalPotential();
  }

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
 * @brief Backtracking algorithm that dynamically adjust step size based on
 * energy evaluation
 * @param f, force object
 * @param dt, initial step size
 * @param rho, discount factor
 * @param time, simulation time
 * @param verbosity, verbosity setting
 * @param totalEnergy_pre, previous energy evaluation
 * @param force, gradient of the energy
 * @param direction, direction, most likely some function of gradient
 * @return
 */
void backtrack(Force &f, const double dt, double rho, double &time,
               const size_t verbosity, const double totalEnergy_pre,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &force,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &direction) {

  // calculate initial energy as reference level
  Eigen::Matrix<double, Eigen::Dynamic, 3> init_position =
      gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);
  double init_time = time;

  // declare variables used in backtracking iterations
  double alpha = dt, totalEnergy, BE, sE, pE, kE, cE, lE, exE;
  size_t count = 0;
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  pos_e += alpha * direction;
  f.update_Vertex_positions();
  std::tie(totalEnergy, BE, sE, pE, kE, cE, lE, exE) = getFreeEnergy(f);

  // while (totalEnergy >
  //        (totalEnergy_pre -
  //         0.02 * alpha * (force.array() * direction.array()).sum())) {
  while (totalEnergy > totalEnergy_pre) {
    if (count > 20) {
      throw std::runtime_error("line search failure!");
    }
    alpha *= rho;
    pos_e = init_position + alpha * direction;
    f.update_Vertex_positions();
    std::tie(totalEnergy, BE, sE, pE, kE, cE, lE, exE) = getFreeEnergy(f);
    // std::cout << "energy pre:" << totalEnergy_pre << std::endl;
    // std::cout << "energy: " << totalEnergy << std::endl;
    count++;
  }

  if (alpha != dt && verbosity > 1) {
    std::cout << "alpha: " << dt << " -> " << alpha << std::endl;
  }
  time = init_time + alpha;
}

/**
 * @brief Save data to richData
 * @param f, force object
 * @param physcialPressure, physical pressre eigen matrix
 * @param verbosity, verbosity setting
 * @return
 */
void saveRichData(
    const Force &f,
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
 * @param energy, components of energy - totalEnergy, BE, sE, pE, kE, cE, lE,
 * exE
 * @param verbosity, verbosity setting
 * @return
 */
void saveNetcdfData(
    const Force &f, size_t &frame, const double &time, TrajFile &fd,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
    const std::tuple<double, double, double, double, double, double, double,
                     double>
        energy,
    const size_t &verbosity) {
  double totalEnergy, BE, sE, pE, kE, cE, lE, exE;
  std::tie(totalEnergy, BE, sE, pE, kE, cE, lE, exE) = energy;
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

  frame = fd.getNextFrameIndex();
  fd.writeTime(frame, time);
  fd.writeCoords(frame, EigenMap<double, 3>(f.vpg.inputVertexPositions));
  fd.writeVelocity(frame, EigenMap<double, 3>(f.vel));
  fd.writeAngles(frame, f.vpg.cornerAngles.raw());

  fd.writeMeanCurvature(frame, H.raw());
  fd.writeSponCurvature(frame, H0.raw());
  fd.writeH_H0_diff(
      frame,
      ((H.raw() - H0.raw()).array() * (H.raw() - H0.raw()).array()).matrix());
  fd.writeExternalPressure(frame, f_ext.raw());
  fd.writePhysicalPressure(frame, fn.raw());
  fd.writeCapillaryPressure(frame, ft.raw());
  fd.writeBendingPressure(frame, fb.raw());
  fd.writeLinePressure(frame, fl.raw());
  fd.writeBendEnergy(frame, BE);
  fd.writeSurfEnergy(frame, sE);
  fd.writePressEnergy(frame, pE);
  fd.writeKineEnergy(frame, kE);
  fd.writeChemEnergy(frame, cE);
  fd.writeLineEnergy(frame, lE);
  fd.writeTotalEnergy(frame, totalEnergy);
}
#endif

void conjugateGradient(Force &f, double dt, double total_time, double tolerance,
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
  Eigen::Matrix<double, Eigen::Dynamic, 3> regularizationForce,
      physicalPressure, DPDForce, direction;

  double totalEnergy, sE, pE, kE, cE, lE, exE,
      oldL2ErrorNorm = 1e6, L2ErrorNorm, dL2ErrorNorm, oldBE = 0.0, BE, dBE,
      dArea, dVolume, dFace, currentNormSq, pastNormSq, time = init_time;

  size_t nMollify = size_t(tMollify / tSave), frame = 0,
         nSave = size_t(tSave / dt);

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
  for (int i = 0; i <= (total_time - init_time) / dt; i++) {

    // compute summerized forces
    getForces(f, physicalPressure, DPDForce, regularizationForce);
    vel_e = physicalPressure + DPDForce + regularizationForce;

    // compute the free energy of the system
    std::tie(totalEnergy, BE, sE, pE, kE, cE, lE, exE) = getFreeEnergy(f);

    // measure the error norm, exit if smaller than tolerance
    L2ErrorNorm = getL2ErrorNorm(physicalPressure);
    if ((i == int((total_time - init_time) / dt)) ||
        (L2ErrorNorm < tolerance)) {
      if (dArea < tolerance && dVolume < tolerance) {
        break;
        if (verbosity > 0) {
          std::cout << "\n"
                    << "Simulation finished, and data saved to " + outputDir
                    << std::endl;
        }
      }else{
        std::cout << "[lambdaSG, lambdaV] = [" << f.P.lambdaSG << ", " << f.P.lambdaV << "]";
        f.P.lambdaSG += f.P.Ksg * (f.surfaceArea - f.targetSurfaceArea) / f.targetSurfaceArea;
        f.P.lambdaV += f.P.Kv * (f.volume - f.refVolume * f.P.Vt) / (f.refVolume * f.P.Vt);
        std::cout << " -> [" << f.P.lambdaSG << ", " << f.P.lambdaV << "]" << std::endl;
      }
    }

    // Save files every nSave iteration and print some info
    if ((i % nSave == 0) || (i == int((total_time - init_time) / dt))) {

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

      // print in-progress information in the console
      if (verbosity > 2) {
        char buffer[50];
        sprintf(buffer, "/t=%d", int(i * dt * 100));
        f.richData.write(outputDir + buffer + ".ply");
        getStatusLog(outputDir + buffer + ".txt", f, dt, i * dt, frame, dArea,
                     dVolume, dBE, dFace, BE, sE, pE, kE, cE, lE, totalEnergy,
                     L2ErrorNorm, f.isTuftedLaplacian, f.isProtein,
                     f.isVertexShift, inputMesh);
      }
      if (verbosity > 1) {
        if (f.P.Ksg != 0 && !f.mesh.hasBoundary()) {
          dArea = abs(f.surfaceArea / f.targetSurfaceArea - 1);
        } else {
          dArea = 0.0;
        }
        if (f.P.Kv != 0 && !f.mesh.hasBoundary()) {
          dVolume = abs(f.volume / f.refVolume / f.P.Vt - 1);
        } else {
          dVolume = 0.0;
        }
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
                  << "\n"
                  << "Increase force spring constant Kf to " << f.P.Kf << "\n";
      }
    }

    // time stepping on vertex position
    if (i % 20 == 0) {
      // pos_e += vel_e * dt;
      // time += dt;
      backtrack(f, dt, 0.5, time, verbosity, totalEnergy, vel_e, vel_e);
      pastNormSq = vel_e.squaredNorm();
      direction = vel_e;
    } else {
      currentNormSq = vel_e.squaredNorm();
      direction = currentNormSq / pastNormSq * direction + vel_e;
      pastNormSq = currentNormSq;
      // pos_e += direction * dt;
      // time += dt;
      backtrack(f, dt, 0.5, time, verbosity, totalEnergy, vel_e, direction);
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
}
} // namespace integration
} // namespace ddgsolver
