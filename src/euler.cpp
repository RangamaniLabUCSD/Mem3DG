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

void backtrack(Force &f, const double dt, double rho, double &time,
               const size_t verbosity, const double totalEnergy_pre,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &force,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &direction);

void euler(Force &f, double dt, double total_time, double tolerance,
           double closeZone, double increment, double maxKv, double maxKsg,
           double tSave, double tMollify, const size_t verbosity,
           std::string inputMesh, std::string outputDir, double init_time,
           double errorJumpLim) {

  // print out a .txt file listing all parameters used
  if (verbosity > 2) {
    getParameterLog(f, dt, total_time, tolerance, tSave, inputMesh, outputDir);
  }

  Eigen::Matrix<double, Eigen::Dynamic, 3> regularizationForce_e;
  regularizationForce_e.resize(f.mesh.nVertices(), 3);
  regularizationForce_e.setZero();

  Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure;
  // , numericalPressure;
  // numericalPressure.resize(f.mesh.nVertices(), 3);
  // numericalPressure.setZero();

  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  // Eigen::Matrix<double, Eigen::Dynamic, 3> pastPosition;
  // pastPosition.resize(f.mesh.nVertices(), 3);
  // pastPosition = pos_e;

  // Eigen::Matrix<double, Eigen::Dynamic, 3> nextPosition;
  // nextPosition.resize(f.mesh.nVertices(), 3);

  // const double hdt = 0.5 * dt, hdt2 = hdt * dt;

  double time, totalEnergy, sE, pE, kE, cE, lE, exE, oldL2ErrorNorm = 1e6, L2ErrorNorm,
                                          dL2ErrorNorm, oldBE = 0.0, BE, dBE,
                                          dArea, dVolume, dFace;
  // double dRef;

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
    time = i * dt + init_time;

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

    vel_e = physicalPressure + regularizationForce_e;

    // periodically save the geometric files, print some info
    if ((i % nSave == 0) || (i == int((total_time - init_time) / dt))) {

      gcs::VertexData<double> H(f.mesh);
      H.fromVector(f.H);
      gcs::VertexData<double> H0(f.mesh);
      H0.fromVector(f.H0);
      gcs::VertexData<double> fn(f.mesh);
      fn.fromVector(rowwiseDotProduct(
          physicalPressure, gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      gcs::VertexData<double> f_ext(f.mesh);
      f_ext.fromVector(
          rowwiseDotProduct(gc::EigenMap<double, 3>(f.externalPressure),
                            gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      gcs::VertexData<double> fb(f.mesh);
      fb.fromVector(
          rowwiseDotProduct(EigenMap<double, 3>(f.bendingPressure),
                            gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      gcs::VertexData<double> fl(f.mesh);
      fl.fromVector(
          rowwiseDotProduct(EigenMap<double, 3>(f.lineTensionPressure),
                            gc::EigenMap<double, 3>(f.vpg.vertexNormals)));
      gcs::VertexData<double> ft(f.mesh);
      ft.fromVector(
          (rowwiseDotProduct(EigenMap<double, 3>(f.capillaryPressure),
                             gc::EigenMap<double, 3>(f.vpg.vertexNormals))
               .array() /
           f.H.array() / 2)
              .matrix());

      std::tie(totalEnergy, BE, sE, pE, kE, cE, lE, exE) = getFreeEnergy(f);

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

      L2ErrorNorm = getL2ErrorNorm(rowwiseScaling(f.mask.cast<double>(), physicalPressure));

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
    } // periodically save the geometric files, print some info

    // integration
    backtrack(f, dt, 0.5, time, verbosity, totalEnergy, vel_e, vel_e);
    //pos_e += vel_e * dt;

    if (f.isProtein) {
      f.proteinDensity.raw() += -f.P.Bc * f.chemicalPotential.raw() * dt;
    }

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
    
  } // integration
} // euler

} // namespace integration
} // namespace ddgsolver
