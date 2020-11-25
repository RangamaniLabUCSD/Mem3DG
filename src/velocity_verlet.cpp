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

void velocityVerlet(Force &f, double dt, double total_time, double tolerance,
                    double closeZone, double increment, double maxKv,
                    double maxKsg, double tSave, double tMollify,
                    const size_t verbosity, std::string inputMesh,
                    std::string outputDir, double init_time,
                    double errorJumpLim) {

  // print out a .txt file listing all parameters used
  if (verbosity > 2) {
    getParameterLog(f, dt, total_time, tolerance, tSave, inputMesh, outputDir);
  }

  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure;
  Eigen::Matrix<double, Eigen::Dynamic, 3> newTotalPressure;
  totalPressure.resize(f.mesh.nVertices(), 3);
  totalPressure.setZero();

  Eigen::Matrix<double, Eigen::Dynamic, 3>  regularizationForce_e;
  regularizationForce_e.resize(f.mesh.nVertices(), 3);
  regularizationForce_e.setZero();

  int nSave = int(tSave / dt);

  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  const double hdt = 0.5 * dt;
  const double hdt2 = hdt * dt;

  Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure;
  Eigen::Matrix<double, Eigen::Dynamic, 3> numericalPressure;

  bool exitFlag = false;

  double totalEnergy;
  double sE;
  double pE;
  double kE;
  double cE;
  double lE;

  double oldL2ErrorNorm = 1e6;
  double L2ErrorNorm;
  double dL2ErrorNorm;

  double oldBE = 0.0;
  double BE;
  double dBE;
  double dArea;
  double dVolume;
  double dFace;
  // double dRef;

  size_t nMollify = size_t(tMollify / tSave);

  std::size_t frame = 0;
#ifdef MEM3DG_WITH_NETCDF
  TrajFile fd;
  if (verbosity > 0) {
    fd.createNewFile(outputDir + "/traj.nc", f.mesh, f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(f.mask.cast<int>());
  }
#endif

  for (int i = 0; i <= (total_time - init_time) / dt; i++) {
    if (f.mesh.hasBoundary()) {
      f.getTubeForces();
    } else {
      f.getVesicleForces();
    }
    f.getDPDForces();
    f.getExternalForces();

    if (f.isProtein) {
      f.getChemicalPotential();
    }

    physicalPressure = gc::EigenMap<double, 3>(f.bendingPressure) +
                       gc::EigenMap<double, 3>(f.capillaryPressure) +
                       gc::EigenMap<double, 3>(f.insidePressure) +
                       gc::EigenMap<double, 3>(f.externalPressure) +
                       gc::EigenMap<double, 3>(f.lineTensionPressure);
    numericalPressure = f.M_inv * (EigenMap<double, 3>(f.dampingForce) +
                                   gc::EigenMap<double, 3>(f.stochasticForce));
    if (!f.mesh.hasBoundary()) {
      removeTranslation(physicalPressure);
      removeRotation(EigenMap<double, 3>(f.vpg.inputVertexPositions),
                     physicalPressure);

      removeTranslation(numericalPressure);
      removeRotation(EigenMap<double, 3>(f.vpg.inputVertexPositions),
                     numericalPressure);
    }

    newTotalPressure = rowwiseScaling(f.mask.cast<double>(),
                                      physicalPressure + numericalPressure);

    // periodically save the geometric files, print some info, compare and
    // adjust
    if ((i % nSave == 0) || (i == int(total_time / dt))) {

      // 1. save
      gcs::VertexData<double> H(f.mesh);
      H.fromVector(f.H);
      f.richData.addVertexProperty("mean_curvature", H);

      gcs::VertexData<double> H0(f.mesh);
      H0.fromVector(f.H0);
      f.richData.addVertexProperty("spon_curvature", H0);

      gcs::VertexData<double> f_ext(f.mesh);
      f_ext.fromVector(f.externalPressureMagnitude);
      f.richData.addVertexProperty("external_pressure", f_ext);

      gcs::VertexData<double> fn(f.mesh);
      fn.fromVector(rowwiseDotProduct(
          physicalPressure, gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      f.richData.addVertexProperty("physical_pressure", fn);

      gcs::VertexData<double> ft(f.mesh);
      ft.fromVector(
          (rowwiseDotProduct(EigenMap<double, 3>(f.capillaryPressure),
                             gc::EigenMap<double, 3>(f.vpg.vertexNormals))
               .array() /
           f.H.array() / 2)
              .matrix());
      f.richData.addVertexProperty("capillary_pressure", ft);

      gcs::VertexData<double> fb(f.mesh);
      fb.fromVector(
          rowwiseDotProduct(EigenMap<double, 3>(f.bendingPressure),
                            gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      f.richData.addVertexProperty("bending_pressure", fb);

      gcs::VertexData<double> fl(f.mesh);
      fl.fromVector(
          rowwiseDotProduct(EigenMap<double, 3>(f.lineTensionPressure),
                            gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      f.richData.addVertexProperty("line_tension_pressure", fl);

      std::tie(totalEnergy, BE, sE, pE, kE, cE, lE) = getFreeEnergy(f);
#ifdef MEM3DG_WITH_NETCDF
      if (verbosity > 0) {
        frame = fd.getNextFrameIndex();
        fd.writeTime(frame, i * dt + init_time);
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

      L2ErrorNorm = getL2ErrorNorm(
          f.M, rowwiseScaling(f.mask.cast<double>(), physicalPressure));
      dL2ErrorNorm = (L2ErrorNorm - oldL2ErrorNorm) / oldL2ErrorNorm;

      if (f.P.Kb != 0) {
        dBE = abs(BE - oldBE) / (BE);
      } else {
        dBE = 0.0;
      }

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

      if (f.P.Ksl != 0 && !f.mesh.hasBoundary()) {
        dFace = ((f.vpg.faceAreas.raw() - f.targetFaceAreas.raw()).array() /
                 f.targetFaceAreas.raw().array())
                    .abs()
                    .sum() /
                f.mesh.nFaces();
      } else {
        dFace = 0.0;
      }

      if (verbosity > 2) {
        char buffer[50];
        sprintf(buffer, "/t=%d", int(i * dt * 100));
        f.richData.write(outputDir + buffer + ".ply");
        getStatusLog(outputDir + buffer + ".txt", f, dt, i * dt, frame, dArea,
                     dVolume, dBE, dFace, BE, sE, pE, kE, cE, lE, totalEnergy,
                     L2ErrorNorm, f.isTuftedLaplacian, f.isProtein,
                     f.isVertexShift, inputMesh);
        // getEnergyLog(i * dt, BE, sE, pE, kE, cE, totalEnergy, outputDir);
      }

      // 2. print
      if (verbosity > 1) {
        std::cout << "\n"
                  << "Time: " << i * dt + init_time << "\n"
                  << "Frame: " << frame << "\n"
                  << "dArea: " << dArea << "\n"
                  << "dVolume:  " << dVolume << "\n"
                  << "dBE: " << dBE << "\n"
                  << "dL2ErrorNorm:   " << dL2ErrorNorm << "\n"
                  << "Bending energy: " << BE << "\n"
                  << "Line energy: " << lE << "\n"
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

      /* for optimization purpose
      // 3.1.1 compare and adjust (in the case of vesicle simulation)
      if ((dVolume < closeZone * tolerance) && (!f.mesh.hasBoundary()) &&
          (dArea < closeZone * tolerance) && (dBE < closeZone * tolerance)) {
        dRef = std::max({dVolume, dArea, dFace});
        if (dRef*increment != 0) {
          f.P.kt *= 1 - dBE / dRef * increment;
          f.P.Kv = std::min(f.P.Kv * (1 + dVolume / dRef * increment), maxKv);
          f.P.Ksg = std::min(f.P.Ksg * (1 + dArea / dRef * increment), maxKsg);
        }

        std::cout << "Within the close zone below " << closeZone
                  << " times tolerance(" << tolerance << "): "
                  << "\n"
                  << "Increase global area penalty Ksg to " << f.P.Ksg << "\n"
                  << "Increase volume penalty Kv to " << f.P.Kv << "\n"
                  << "Decrese randomness kT to " << f.P.kt << "\n";
      }

      // 3.1.2 increase the force spring constant
      f.P.Kf *= 1 + increment;

      // 3.2 compare and exit
      if (((dVolume < tolerance) && (dArea < tolerance) && (dBE < tolerance) &&
      (dL2ErrorNorm > 0))||(exitFlag == true)) { f.P.Kv *= (1 - dVolume / dRef *
      increment); f.P.Ksg *= (1 - dArea / dRef * increment); exitFlag = true;
        f.P.kt = 0.0;
        increment = 0;
        if (nMollify > 0) {

          std::cout << "\n"
                    << nMollify << " mollification(s) left" << std::endl;
          nMollify -= 1;
        } else {
          std::cout << "\n"
                    << "Converged! Saved to " + outputDir << std::endl;
          f.richData.write(outputDir + "final.ply");
          getStatusLog(outputDir + "summary.txt", f, dt, i * dt, frame, dArea,
                       dVolume, dBE, dFace, BE, sE, pE, cE, totalEnergy,
                       L2ErrorNorm, f.isTuftedLaplacian, f.isProtein,
      f.isVertexShift, inputMesh); break;
        }
      }
      */

      // 3.3 fail and exit
      if (abs(dL2ErrorNorm) > errorJumpLim) {
        if (verbosity > 0) {
          std::cout << "Error Norm changes rapidly. Save data and quit."
                    << std::endl;
        }
        break;
      }

      oldL2ErrorNorm = L2ErrorNorm;
      oldBE = BE;

    } // periodically save the geometric files, print some info, compare and
      // adjust

    // integration
    pos_e += vel_e * dt + hdt2 * totalPressure;
    vel_e += (totalPressure + newTotalPressure) * hdt;
    totalPressure = newTotalPressure;

    if (f.isProtein) {
      f.proteinDensity.raw() += -f.P.Bc * f.chemicalPotential.raw() * dt;
    }

    // Regularize the vetex position geometry if needed
    regularizationForce_e = rowwiseScaling(f.mask.cast<double>(),
                                      gc::EigenMap<double, 3>(f.regularizationForce));
    pos_e +=  regularizationForce_e * dt;

    if (f.isVertexShift) {
      vertexShift(f.mesh, f.vpg, f.mask);
    }

    // if (f.vpg.cornerAngles.raw().minCoeff() < (M_PI / 6)) {
    //   f.isTuftedLaplacian = true;
    // } else {
    //   f.isTuftedLaplacian = false;
    // }

    // recompute cached values
    f.update_Vertex_positions();

    /* for optimization purpose
    // 3.3.A fail and exit
    if (i == int(total_time / dt)) {
      std::cout
          << "\n"
          << "Fail to converge in given time and Exit. Past data saved to " +
                 outputDir
          << std::endl;
      getStatusLog(outputDir + "failure_report.txt", f, dt, i * dt, frame,
                   dArea, dVolume, dBE, dFace, BE, sE, pE, cE, totalEnergy,
                   L2ErrorNorm, f.isTuftedLaplacian, f.isProtein,
    f.isVertexShift, inputMesh);
    }
    */

    // 3.3.B finish and exit
    if (verbosity > 0) {
      if (i == int((total_time - init_time) / dt)) {
        std::cout << "\n"
                  << "Simulation finished, and data saved to " + outputDir
                  << std::endl;
        getStatusLog(outputDir + "/final_report.txt", f, dt, i * dt, frame,
                     dArea, dVolume, dBE, dFace, BE, sE, pE, kE, cE, lE,
                     totalEnergy, L2ErrorNorm, f.isTuftedLaplacian, f.isProtein,
                     f.isVertexShift, inputMesh);
      }
    }
  } // periodic save, print and adjust
}

} // namespace integration
} // namespace ddgsolver
