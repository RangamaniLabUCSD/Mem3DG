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
 * @param EXIT, exit flag for integration loop
 * @param verbosity, verbosity setting
 * @param totalEnergy_pre, previous energy evaluation
 * @param force, gradient of the energy
 * @param direction, direction, most likely some function of gradient
 * @return
 */
void backtrack(System &f, const double dt, double rho, double &time, bool &EXIT,
               const size_t verbosity, const double totalEnergy_pre,
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

  while (f.E.totalE >
         (totalEnergy_pre -
          0 * alpha * (force.array() * direction.array()).sum())) {
    // while (f.E.totalE > totalEnergy_pre) {
    if (count > 50) {
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
    // std::cout << "energy pre:" << totalEnergy_pre << std::endl;
    // std::cout << "energy: " << f.E.totalE << std::endl;
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
    const std::tuple<double, double, double, double, double, double, double,
                     double>
        energy,
    const size_t &verbosity) {

  Eigen::Matrix<double, Eigen::Dynamic, 1> fn, f_ext, fb, fl, ft;
  double totalE, BE, sE, pE, kE, cE, lE, exE;
  std::tie(totalE, BE, sE, pE, kE, cE, lE, exE) = energy;

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

  fd.writeMeanCurvature(frame, f.H);
  fd.writeSponCurvature(frame, f.H0);
  fd.writeH_H0_diff(frame,
                    ((f.H - f.H0).array() * (f.H - f.H0).array()).matrix());
  fd.writeExternalPressure(frame, f_ext);
  fd.writePhysicalPressure(frame, fn);
  fd.writeCapillaryPressure(frame, ft);
  fd.writeBendingPressure(frame, fb);
  fd.writeLinePressure(frame, fl);
  fd.writeBendEnergy(frame, BE);
  fd.writeSurfEnergy(frame, sE);
  fd.writePressEnergy(frame, pE);
  fd.writeKineEnergy(frame, kE);
  fd.writeChemEnergy(frame, cE);
  fd.writeLineEnergy(frame, lE);
  fd.writeTotalEnergy(frame, totalE);
}
#endif

void conjugateGradient(System &f, double dt, double total_time,
                       double tolerance, double closeZone, double increment,
                       double maxKv, double maxKsg, double tSave,
                       double tMollify, const size_t verbosity,
                       std::string inputMesh, std::string outputDir,
                       double init_time, double errorJumpLim) {

  // print out a txt file listing all parameters used
  if (verbosity > 2) {
    getParameterLog(f, dt, total_time, tolerance, tSave, inputMesh, outputDir);
  }

  // initialize variables used in time integration
  Eigen::Matrix<double, Eigen::Dynamic, 3> regularizationForce,
      physicalPressure, DPDForce, direction;
  double dArea, dVolume, currentNormSq, pastNormSq,
      time = init_time;
  size_t frame = 0;
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
    vel_e = physicalPressure + DPDForce + regularizationForce;

    // measure the error norm and constraints, exit if smaller than tolerance
    dArea = (f.P.Ksg != 0 && !f.mesh.hasBoundary())
                ? abs(f.surfaceArea / f.targetSurfaceArea - 1)
                : 0.0;
    dVolume = (f.P.Kv != 0 && !f.mesh.hasBoundary())
                  ? abs(f.volume / f.refVolume / f.P.Vt - 1)
                  : 0.0;
    f.getL2ErrorNorm(physicalPressure);
    if (f.L2ErrorNorm < tolerance) {
      if (dArea < tolerance && dVolume < tolerance) {
        EXIT = true;
      } else {
        std::cout << "\n[lambdaSG, lambdaV] = [" << f.P.lambdaSG << ", "
                  << f.P.lambdaV << "]";
        f.P.lambdaSG += f.P.Ksg * (f.surfaceArea - f.targetSurfaceArea) /
                        f.targetSurfaceArea;
        f.P.lambdaV +=
            f.P.Kv * (f.volume - f.refVolume * f.P.Vt) / (f.refVolume * f.P.Vt);
        std::cout << " -> [" << f.P.lambdaSG << ", " << f.P.lambdaV << "]"
                  << std::endl;
      }
    }

    // compute the free energy of the system
    f.getFreeEnergy();

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
                       std::tie(f.E.totalE, f.E.BE, f.E.sE, f.E.pE, f.E.kE, f.E.cE, f.E.lE, f.E.exE),
                       verbosity);
      }
#endif

      // print in-progress information in the console
      if (verbosity > 2) {
        char buffer[50];
        sprintf(buffer, "/t=%d", int(time * 100));
        f.richData.write(outputDir + buffer + ".ply");
      }
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
        std::cout << "\n"
                  << "Simulation finished, and data saved to " + outputDir
                  << std::endl;
      }
      break;
    }

    // time stepping on vertex position
    size_t countCG = 0;
    if (countCG % 20 == 0) {
      pastNormSq = vel_e.squaredNorm();
      direction = vel_e;
      countCG = 0;
    } else {
      currentNormSq = vel_e.squaredNorm();
      direction = currentNormSq / pastNormSq * direction + vel_e;
      pastNormSq = currentNormSq;
      countCG++;
    }
    // pos_e += direction * dt;
    // time += dt;
    backtrack(f, dt, 0.5, time, EXIT, verbosity, f.E.totalE, vel_e, direction);
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
