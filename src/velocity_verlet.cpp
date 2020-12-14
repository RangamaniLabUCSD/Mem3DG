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

void getForces(Force &f,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &DPDForce,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &regularizationForce);

void velocityVerlet(Force &f, double dt, double total_time, double tolerance,
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
  double totalEnergy, sE, pE, kE, cE, lE, exE,
      oldL2ErrorNorm = 1e6, L2ErrorNorm, dL2ErrorNorm, oldBE = 0.0, BE, dBE,
      dArea, dVolume, dFace; // double dRef;

  size_t nMollify = size_t(tMollify / tSave), frame = 0,
         nSave = size_t(tSave / dt);

  // map the raw eigen datatype for computation
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg.inputVertexPositions);

#ifdef MEM3DG_WITH_NETCDF
  TrajFile fd;
  if (verbosity > 0) {
    fd.createNewFile(outputDir + "/traj.nc", f.mesh, f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(f.mask.cast<int>());
  }
#endif

  for (int i = 0; i <= (total_time - init_time) / dt; i++) {

    // compute summerized forces
    getForces(f, physicalPressure, DPDForce, regularizationForce);
    totalPressure.resize(f.mesh.nVertices(), 3);
    totalPressure.setZero();
    newTotalPressure = physicalPressure + DPDForce;

    // measure the error norm, exit if smaller than tolerance or reach time
    // limit
    L2ErrorNorm = getL2ErrorNorm(physicalPressure);
    if ((i == int((total_time - init_time) / dt)) || (L2ErrorNorm < 1e-3)) {
      break;
      if (verbosity > 0) {
        std::cout << "\n"
                  << "Simulation finished, and data saved to " + outputDir
                  << std::endl;
      }
    }

    // Save files every nSave iteration and print some info
    if ((i % nSave == 0) || (i == int((total_time - init_time) / dt))) {

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
      ft.fromVector(
          (rowwiseDotProduct(EigenMap<double, 3>(f.capillaryPressure),
                             gc::EigenMap<double, 3>(f.vpg.vertexNormals))
               .array() /
           f.H.array() / 2)
              .matrix());

      std::tie(totalEnergy, BE, sE, pE, kE, cE, lE, exE) = getFreeEnergy(f);

      // save variable to richData
      if (verbosity > 2) {
        f.richData.addVertexProperty("mean_curvature", H);
        f.richData.addVertexProperty("spon_curvature", H0);
        f.richData.addVertexProperty("external_pressure", f_ext);
        f.richData.addVertexProperty("physical_pressure", fn);
        f.richData.addVertexProperty("capillary_pressure", ft);
        f.richData.addVertexProperty("bending_pressure", fb);
        f.richData.addVertexProperty("line_tension_pressure", fl);
      }

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

      // print in-progress information in the console
      dL2ErrorNorm = (L2ErrorNorm - oldL2ErrorNorm) / oldL2ErrorNorm;
      if (abs(dL2ErrorNorm) > errorJumpLim) {
        if (verbosity > 0) {
          std::cout << "Error Norm changes rapidly. Save data and quit."
                    << std::endl;
        }
        break;
      }
      oldL2ErrorNorm = L2ErrorNorm;
      oldBE = BE;

      if (verbosity > 1) {
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

      if (verbosity > 2) {
        char buffer[50];
        sprintf(buffer, "/t=%d", int(i * dt * 100));
        f.richData.write(outputDir + buffer + ".ply");
        getStatusLog(outputDir + buffer + ".txt", f, dt, i * dt, frame, dArea,
                     dVolume, dBE, dFace, BE, sE, pE, kE, cE, lE, totalEnergy,
                     L2ErrorNorm, f.isTuftedLaplacian, f.isProtein,
                     f.isVertexShift, inputMesh);
      }
    }

    // time stepping on vertex position
    pos_e += vel_e * dt + hdt2 * totalPressure;
    pos_e += regularizationForce * dt;
    vel_e += (totalPressure + newTotalPressure) * hdt;
    totalPressure = newTotalPressure;
    if (f.isVertexShift) {
      vertexShift(f.mesh, f.vpg, f.mask);
    }

    // time stepping on protein density
    if (f.isProtein) {
      f.proteinDensity.raw() += -f.P.Bc * f.chemicalPotential.raw() * dt;
    }

    // recompute cached values
    f.update_Vertex_positions();
  }
}

} // namespace integration
} // namespace ddgsolver
