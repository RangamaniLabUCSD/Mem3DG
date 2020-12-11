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

Eigen::Matrix<double, Eigen::Dynamic, 3>
getForces(Force &f, Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
          Eigen::Matrix<double, Eigen::Dynamic, 3> &regularizationForce_e) {
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

  regularizationForce_e = rowwiseScaling(
      f.mask.cast<double>(), gc::EigenMap<double, 3>(f.regularizationForce));

  // numericalPressure = f.M_inv * (EigenMap<double, 3>(f.dampingForce) +
  //                                gc::EigenMap<double,
  //                                3>(f.stochasticForce));

  if (!f.mesh.hasBoundary()) {
    removeTranslation(physicalPressure);
    removeRotation(EigenMap<double, 3>(f.vpg.inputVertexPositions),
                   physicalPressure);

    // removeTranslation(numericalPressure);
    // removeRotation(EigenMap<double, 3>(f.vpg.inputVertexPositions),
    //                numericalPressure);
  }

  return physicalPressure + regularizationForce_e;
}

void backtrack(Force &f, const double dt, double rho, double &time,
               const size_t verbosity, const double totalEnergy_pre,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &force,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &direction) {

  // calculate initial energy as reference level
  Eigen::Matrix<double, Eigen::Dynamic, 3> init_position =
      gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);
  double init_time = time;

  // declare variables used in backtracking iterations
  double alpha = dt, totalEnergy, BE, sE, pE, kE, cE, lE;
  size_t count = 0;
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  pos_e += alpha * direction;
  f.update_Vertex_positions();
  std::tie(totalEnergy, BE, sE, pE, kE, cE, lE) = getFreeEnergy(f);

  // while (totalEnergy >
  //        (totalEnergy_pre -
  //         0.05 * alpha * (force.array() * direction.array()).sum())) {
  while (totalEnergy > totalEnergy_pre) {
    if (count > 20) {
      throw std::runtime_error("line search failure!");
    }
    alpha *= rho;
    pos_e = init_position + alpha * direction;
    f.update_Vertex_positions();
    std::tie(totalEnergy, BE, sE, pE, kE, cE, lE) = getFreeEnergy(f);
    std::cout << totalEnergy - totalEnergy_pre << std::endl;
    // std::cout << totalEnergy_pre << std::endl;
    //std::cout << totalEnergy << std::endl;
    count++;
  }

  if (alpha != dt && verbosity > 1){
    std::cout << "dt = " << dt << " -> " << alpha << std::endl;
  }
  time = init_time + alpha;
}

void conjugateGradient(Force &f, double dt, double total_time, double tolerance,
                       double closeZone, double increment, double maxKv,
                       double maxKsg, double tSave, double tMollify,
                       const size_t verbosity, std::string inputMesh,
                       std::string outputDir, double init_time,
                       double errorJumpLim) {

  // print out a .txt file listing all parameters used
  if (verbosity > 2) {
    getParameterLog(f, dt, total_time, tolerance, tSave, inputMesh, outputDir);
  }

  Eigen::Matrix<double, Eigen::Dynamic, 3> regularizationForce_e;
  regularizationForce_e.resize(f.mesh.nVertices(), 3);
  regularizationForce_e.setZero();

  Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure, direction;

  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  double totalEnergy, sE, pE, kE, cE, lE,
      oldL2ErrorNorm = 1e6, L2ErrorNorm = 1e6, dL2ErrorNorm, oldBE = 0.0, BE,
      dBE, dArea, dVolume, dFace, currentNorm2, pastNorm2, time = init_time;

  size_t nMollify = size_t(tMollify / tSave), frame = 0,
         nSave = size_t(tSave / dt);

#ifdef MEM3DG_WITH_NETCDF
  TrajFile fd;
  if (verbosity > 0) {
    fd.createNewFile(outputDir + "/traj.nc", f.mesh, f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(f.mask.cast<int>());
  }
#endif

  for (int i = 0; i <= (total_time - init_time) / dt; i++) {

    vel_e = getForces(f, physicalPressure, regularizationForce_e);
    L2ErrorNorm = getL2ErrorNorm(rowwiseScaling(f.mask.cast<double>(), physicalPressure));

    if ((i == int((total_time - init_time) / dt)) || (L2ErrorNorm < 1e-3)) {
      break;
      if (verbosity > 0) {
        std::cout << "\n"
                  << "Simulation finished, and data saved to " + outputDir
                  << std::endl;
      }
    }

    // 1. save
    // periodically save the geometric files, print some info
    if ((i % nSave == 0) || (i == int((total_time - init_time) / dt))) {

      gcs::VertexData<double> H(f.mesh);
      H.fromVector(f.H);
      gcs::VertexData<double> H0(f.mesh);
      H0.fromVector(f.H0);
      gcs::VertexData<double> f_ext(f.mesh);
      f_ext.fromVector(f.externalPressureMagnitude);
      gcs::VertexData<double> fn(f.mesh);
      fn.fromVector(rowwiseDotProduct(
          physicalPressure, gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      gcs::VertexData<double> ft(f.mesh);
      ft.fromVector(
          (rowwiseDotProduct(EigenMap<double, 3>(f.capillaryPressure),
                             gc::EigenMap<double, 3>(f.vpg.vertexNormals))
               .array() /
           f.H.array() / 2)
              .matrix());
      gcs::VertexData<double> fb(f.mesh);
      fb.fromVector(
          rowwiseDotProduct(EigenMap<double, 3>(f.bendingPressure),
                            gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      gcs::VertexData<double> fl(f.mesh);
      fl.fromVector(
          rowwiseDotProduct(EigenMap<double, 3>(f.lineTensionPressure),
                            gc::EigenMap<double, 3>(f.vpg.vertexNormals)));

      std::tie(totalEnergy, BE, sE, pE, kE, cE, lE) = getFreeEnergy(f);

      // 1.1 save to .ply file
      if (verbosity > 2) {
        f.richData.addVertexProperty("mean_curvature", H);
        f.richData.addVertexProperty("spon_curvature", H0);
        f.richData.addVertexProperty("external_pressure", f_ext);
        f.richData.addVertexProperty("physical_pressure", fn);
        f.richData.addVertexProperty("capillary_pressure", ft);
        f.richData.addVertexProperty("bending_pressure", fb);
        f.richData.addVertexProperty("line_tension_pressure", fl);
      }

      // 1.2 save to .nc file
#ifdef MEM3DG_WITH_NETCDF
      if (verbosity > 0) {
        frame = fd.getNextFrameIndex();
        fd.writeTime(frame, time);
        fd.writeCoords(frame, EigenMap<double, 3>(f.vpg.inputVertexPositions));
        fd.writeVelocity(frame, EigenMap<double, 3>(f.vel));
        fd.writeAngles(frame, f.vpg.cornerAngles.raw());

        fd.writeMeanCurvature(frame, H.raw());
        fd.writeSponCurvature(frame, H0.raw());
        fd.writeH_H0_diff(
            frame, ((H.raw() - H0.raw()).array() * (H.raw() - H0.raw()).array())
                       .matrix());
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

      if (verbosity > 2) {
        char buffer[50];
        sprintf(buffer, "/t=%d", int(i * dt * 100));
        f.richData.write(outputDir + buffer + ".ply");
        getStatusLog(outputDir + buffer + ".txt", f, dt, i * dt, frame, dArea,
                     dVolume, dBE, dFace, BE, sE, pE, kE, cE, lE, totalEnergy,
                     L2ErrorNorm, f.isTuftedLaplacian, f.isProtein,
                     f.isVertexShift, inputMesh);
      }

      // 2. print
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

    /// time integration
    if (i % 20 == 0) {
      pos_e += vel_e * dt;
      pastNorm2 = vel_e.squaredNorm();
      direction = vel_e;
      time += dt;
    } else {
      currentNorm2 = vel_e.squaredNorm();
      direction = currentNorm2 / pastNorm2 * direction + vel_e;
      pastNorm2 = currentNorm2;
      // pos_e += direction * dt;
      backtrack(f, dt, 0.5, time, verbosity, totalEnergy, vel_e, direction);
      std::tie(totalEnergy, BE, sE, pE, kE, cE, lE) = getFreeEnergy(f);
    }

    if (f.isProtein) {
      f.proteinDensity.raw() += -f.P.Bc * f.chemicalPotential.raw() * dt;
    }

    if (f.isVertexShift) {
      vertexShift(f.mesh, f.vpg, f.mask);
    }

    // recompute cached values
    f.update_Vertex_positions();
    
  } // integration 
} // euler

} // namespace integration
} // namespace ddgsolver
