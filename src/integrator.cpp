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

#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/solver/util.h"
#include "mem3dg/solver/version.h"

#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace std;

namespace mem3dg {

double Integrator::backtrack(
    double rho, double c1, bool &EXIT, bool &SUCCESS,
    const double potentialEnergy_pre,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &physicalForce,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &direction) {

  // calculate initial energy as reference level
  Eigen::Matrix<double, Eigen::Dynamic, 3> init_position =
      gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  double init_time = f.time;

  // declare variables used in backtracking iterations
  double alpha = dt;
  size_t count = 0;
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);

  pos_e += alpha * direction;
  f.updateVertexPositions();
  f.computeFreeEnergy();
  double projection =
      (physicalForce.array() *
       rowwiseDotProduct(direction, EigenMap<double, 3>(f.vpg->vertexNormals))
           .array())
          .sum();

  while (true) {
    if (projection < 0) {
      std::cout << "\nBacktracking line search: on uphill direction! \n"
                << std::endl;
      EXIT = true;
      SUCCESS = false;
      break;
    }
    // while (f.E.potE > potentialEnergy_pre) {
    if (f.E.potE < (potentialEnergy_pre - c1 * alpha * projection)) {
      break;
    }
    if (alpha < 1e-12) {
      std::cout << "\nline search failure! Simulation stopped. \n" << std::endl;
      EXIT = true;
      SUCCESS = false;
      break;
    }
    alpha *= rho;
    pos_e = init_position + alpha * direction;
    f.updateVertexPositions();
    f.computeFreeEnergy();
    count++;
  }

  if (alpha != dt && verbosity > 3) {
    std::cout << "alpha: " << dt << " -> " << alpha << std::endl;
    std::cout << "L1 norm: " << f.L1ErrorNorm << std::endl;
  }
  f.time = init_time + alpha;

  return alpha;
}

void Integrator::getForces() {
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(f.vpg->vertexNormals);

  f.computePhysicalForces();

  physicalForce =
      (f.mask.raw().cast<double>()).array() *
      (f.bendingForce.raw() + f.capillaryForce.raw() + f.externalForce.raw() +
       f.osmoticForce.raw() + f.lineCapillaryForce.raw())
          .array();

  if ((f.P.gamma != 0) || (f.P.temp != 0)) {
    f.computeDPDForces(dt);
    DPDForce = f.mask.raw().cast<double>().array() *
               rowwiseDotProduct((EigenMap<double, 3>(f.dampingForce) +
                                  EigenMap<double, 3>(f.stochasticForce)),
                                 vertexAngleNormal_e)
                   .array();
  }

  if (!f.mesh->hasBoundary()) {
    removeTranslation(rowwiseScaling(physicalForce, vertexAngleNormal_e));
    removeRotation(EigenMap<double, 3>(f.vpg->inputVertexPositions),
                   rowwiseScaling(physicalForce, vertexAngleNormal_e));
    // removeTranslation(DPDPressure);
    // removeRotation(EigenMap<double, 3>(f.vpg->inputVertexPositions),
    // DPDPressure);
  }
}

void Integrator::pressureConstraintThreshold(bool &EXIT,
                                             const bool isAugmentedLagrangian,
                                             const double dArea,
                                             const double ctol,
                                             double increment) {
  if (f.L1ErrorNorm < tol) {
    if (isAugmentedLagrangian) { // augmented Lagrangian method
      if (dArea < ctol) {        // exit if fulfilled all constraints
        std::cout << "\nL1 error norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[lambdaSG] = [" << f.P.lambdaSG << ", "
                  << "]";
        f.P.lambdaSG +=
            f.P.Ksg * (f.surfaceArea - f.refSurfaceArea) / f.refSurfaceArea;
        std::cout << " -> [" << f.P.lambdaSG << "]" << std::endl;
      }
    } else {              // incremental harmonic penalty method
      if (dArea < ctol) { // exit if fulfilled all constraints
        std::cout << "\nL1 error norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[Ksg] = [" << f.P.Ksg << "]";
        f.P.Ksg *= increment;
        std::cout << " -> [" << f.P.Ksg << "]" << std::endl;
      }
    }
  }
}

void Integrator::reducedVolumeThreshold(bool &EXIT,
                                        const bool isAugmentedLagrangian,
                                        const double dArea,
                                        const double dVolume, const double ctol,
                                        double increment) {
  if (f.L1ErrorNorm < tol) {
    if (isAugmentedLagrangian) {            // augmented Lagrangian method
      if (dArea < ctol && dVolume < ctol) { // exit if fulfilled all constraints
        std::cout << "\nL1 error norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[lambdaSG, lambdaV] = [" << f.P.lambdaSG << ", "
                  << f.P.lambdaV << "]";
        f.P.lambdaSG +=
            f.P.Ksg * (f.surfaceArea - f.refSurfaceArea) / f.refSurfaceArea;
        f.P.lambdaV +=
            f.P.Kv * (f.volume - f.refVolume * f.P.Vt) / (f.refVolume * f.P.Vt);
        std::cout << " -> [" << f.P.lambdaSG << ", " << f.P.lambdaV << "]"
                  << std::endl;
      }
    } else { // incremental harmonic penalty method
      if (dArea < ctol && dVolume < ctol) { // exit if fulfilled all constraints
        std::cout << "\nL1 error norm smaller than tolerance." << std::endl;
        EXIT = true;
      }

      // iterate if not
      if (dArea > ctol) {
        std::cout << "\n[Ksg] = [" << f.P.Ksg << "]";
        f.P.Ksg *= 1.3;
        std::cout << " -> [" << f.P.Ksg << "]" << std::endl;
      }
      if (dVolume > ctol) {
        std::cout << "\n[Kv] = [" << f.P.Kv << "]";
        f.P.Kv *= 1.3;
        std::cout << " -> [" << f.P.Kv << "]" << std::endl;
      }
    }
  }
}

void Integrator::createNetcdfFile() {
  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  if (verbosity > 0) {
    fd.createNewFile(outputDir + trajFileName, *f.mesh, *f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(f.mask.raw().cast<int>());
    if (!f.mesh->hasBoundary()) {
      fd.writeRefVolume(f.refVolume);
      fd.writeRefSurfArea(f.refSurfaceArea);
    }
  }
#endif
}

void Integrator::saveData() {
  // save variable to richData and save ply file
  if ((verbosity > 3 && !f.O.isGrowMesh) || (verbosity > 0 && f.O.isGrowMesh)) {
    char buffer[50];
    sprintf(buffer, "/frame%d.ply", (int)frame);
    saveRichData(buffer);
  }

#ifdef MEM3DG_WITH_NETCDF
  // save variable to netcdf traj file
  if (verbosity > 0) {
    saveNetcdfData();
  }
#endif

  // print in-progress information in the console
  if (verbosity > 1) {
    std::cout << "\n"
              << "t: " << f.time << ", "
              << "n: " << frame << "\n"
              << "dA/Area: " << dArea << "/" << f.surfaceArea << ", "
              << "dVP/Volume: " << dVP << "/" << f.volume << ", "
              << "h: "
              << abs(f.vpg->inputVertexPositions[f.thePoint.nearestVertex()].z)
              << "\n"
              << "E_total: " << f.E.totalE << "\n"
              << "E_pot: " << f.E.potE << "\n"
              << "|e|L1: " << f.L1ErrorNorm << "\n"
              << "H: ["
              << (f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                  f.vpg->vertexMeanCurvatures.raw())
                     .minCoeff()
              << ","
              << (f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                  f.vpg->vertexMeanCurvatures.raw())
                     .maxCoeff()
              << "]"
              << "\n"
              << "K: ["
              << (f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                  f.vpg->vertexGaussianCurvatures.raw())
                     .minCoeff()
              << ","
              << (f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                  f.vpg->vertexGaussianCurvatures.raw())
                     .maxCoeff()
              << "]" << std::endl;
    // << "COM: "
    // << gc::EigenMap<double,
    // 3>(f.vpg->inputVertexPositions).colwise().sum() /
    //         f.vpg->inputVertexPositions.raw().rows()
    // << "\n"
  }
  // break loop if EXIT flag is on
  if (EXIT) {
    if (verbosity > 0) {
      std::cout << "Simulation " << (SUCCESS ? "finished" : "failed")
                << ", and data saved to " + outputDir << std::endl;
      if (verbosity > 2) {
        saveRichData("/out.ply");
      }
    }
  }

  frame++;
}

void Integrator::markFileName(std::string marker_str) {
  std::string dirPath = outputDir;

  const char *marker = marker_str.c_str();

  char *file = new char[trajFileName.size() + 1];
  std::copy(trajFileName.begin(), trajFileName.end(), file);
  file[trajFileName.size()] = '\0';

  char fileMarked[50]{"/"}, oldNC[150]{"/"}, newNC[150]{"/"};

  // sprintf(fileMarked, "/traj_H_%d_VP_%d_failed.nc", int(H * 100),
  //         int(VP * 100));

  // split the extension and file name
  const char *ext = strchr(file, '.');

  // name fileMarked to be the file name
  strncpy(fileMarked, file, ext - file);

  // name fileMarked to be file name + the marker + extension
  strcat(fileMarked, marker);
  strcat(fileMarked, ext);
  fileMarked[ext - file + sizeof(marker) + sizeof(ext)] = '\0';

  // append the directory path and copy to oldNC and newNC
  strcpy(oldNC, dirPath.c_str());
  strcpy(newNC, dirPath.c_str());
  strcat(oldNC, file);
  strcat(newNC, fileMarked);

  // rename file
  rename(oldNC, newNC);
  delete[] file;
}

void Integrator::saveRichData(std::string plyName) {
  gcs::RichSurfaceMeshData richData(*f.mesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*f.vpg);

  // write protein distribution
  if (f.O.isProtein) {
    richData.addVertexProperty("protein_density", f.proteinDensity);
  }

  // write bool
  gcs::VertexData<int> msk(*f.mesh);
  msk.fromVector(f.mask.raw().cast<int>());
  richData.addVertexProperty("mask", msk);
  gcs::VertexData<int> tkr(*f.mesh);
  tkr.fromVector(f.thePointTracker.raw().cast<int>());
  richData.addVertexProperty("the_point", tkr);

  // write geometry
  gcs::VertexData<double> meanCurv(*f.mesh);
  meanCurv.fromVector(f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                      f.vpg->vertexMeanCurvatures.raw());
  richData.addVertexProperty("mean_curvature", meanCurv);
  gcs::VertexData<double> gaussCurv(*f.mesh);
  gaussCurv.fromVector(f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                       f.vpg->vertexGaussianCurvatures.raw());
  richData.addVertexProperty("gauss_curvature", gaussCurv);
  richData.addVertexProperty("spon_curvature", f.H0);

  // write pressures
  gcs::VertexData<double> fn(*f.mesh);
  fn.fromVector(physicalForce);
  richData.addVertexProperty("bending_force", f.bendingForce);
  richData.addVertexProperty("capillary_force", f.capillaryForce);
  richData.addVertexProperty("line_tension_force", f.lineCapillaryForce);
  richData.addVertexProperty("osmotic_force", f.osmoticForce);
  richData.addVertexProperty("external_force", f.externalForce);
  richData.addVertexProperty("physical_force", fn);

  richData.write(outputDir + plyName);
}

#ifdef MEM3DG_WITH_NETCDF
void Integrator::saveNetcdfData() {
  std::size_t idx = fd.getNextFrameIndex();

  // scalar quantities
  // write time
  fd.writeTime(idx, f.time);
  // write geometry
  fd.writeVolume(idx, f.volume);
  fd.writeSurfArea(idx, f.mesh->hasBoundary() ? f.surfaceArea - f.refSurfaceArea
                                              : f.surfaceArea);
  fd.writeHeight(
      idx, abs(f.vpg->inputVertexPositions[f.thePoint.nearestVertex()].z));
  // write energies
  fd.writeBendEnergy(idx, f.E.BE);
  fd.writeSurfEnergy(idx, f.E.sE);
  fd.writePressEnergy(idx, f.E.pE);
  fd.writeKineEnergy(idx, f.E.kE);
  fd.writeChemEnergy(idx, f.E.cE);
  fd.writeLineEnergy(idx, f.E.lE);
  fd.writeTotalEnergy(idx, f.E.totalE);
  // write Norms
  fd.writeL1ErrorNorm(idx, f.L1ErrorNorm);
  fd.writeL1BendNorm(idx, f.computeL1Norm(f.bendingForce.raw()));
  fd.writeL1SurfNorm(idx, f.computeL1Norm(f.capillaryForce.raw()));
  fd.writeL1PressNorm(idx, f.computeL1Norm(f.osmoticForce.raw()));
  fd.writeL1LineNorm(idx, f.computeL1Norm(f.lineCapillaryForce.raw()));

  // vector quantities
  if (!f.O.isGrowMesh) {
    // write velocity
    fd.writeVelocity(idx, EigenMap<double, 3>(f.vel));
    // write protein density distribution
    if (f.O.isProtein) {
      fd.writeProteinDensity(idx, f.proteinDensity.raw());
    }

    // write geometry
    fd.writeCoords(idx, EigenMap<double, 3>(f.vpg->inputVertexPositions));
    fd.writeTopoFrame(idx, getFaceVertexMatrix(*f.mesh));
    fd.writeMeanCurvature(idx, f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                                   f.vpg->vertexMeanCurvatures.raw());
    fd.writeGaussCurvature(idx, f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                                    f.vpg->vertexGaussianCurvatures.raw());
    fd.writeSponCurvature(idx, f.H0.raw());
    // fd.writeAngles(idx, f.vpg.cornerAngles.raw());
    // fd.writeH_H0_diff(idx,
    //                   ((f.H - f.H0).array() * (f.H -
    //                   f.H0).array()).matrix());

    // write pressures
    fd.writeBendingForce(idx, f.bendingForce.raw());
    fd.writeCapillaryForce(idx, f.capillaryForce.raw());
    fd.writeLineForce(idx, f.lineCapillaryForce.raw());
    fd.writeOsmoticForce(idx, f.osmoticForce.raw());
    fd.writeExternalForce(idx, f.externalForce.raw());
    fd.writePhysicalForce(idx, physicalForce);
  }
}
#endif

void Integrator::getParameterLog(std::string inputMesh) {
  ofstream myfile(outputDir + "/parameter.txt");
  if (myfile.is_open()) {
    myfile << "Mem3DG Version: " << MEM3DG_VERSION << "\n";
    myfile << "Input Mesh:     " << inputMesh << "\n";
    myfile << "Physical parameters used: \n";
    myfile << "\n";
    myfile << "Kb(bare):     " << f.P.Kb << "\n"
           << "Kb(coated):   " << f.P.Kbc << "\n"
           << "H0:     " << f.P.H0 << "\n"
           << "Kse:    " << f.P.Kse << "\n"
           << "Ksl:    " << f.P.Ksl << "\n"
           << "Kst:    " << f.P.Kst << "\n"
           << "Ksg:    " << f.P.Ksg << "\n"
           << "Kv:     " << f.P.Kv << "\n"
           << "gamma:  " << f.P.gamma << "\n"
           << "Vt:     " << f.P.Vt << "\n"
           << "kt:     " << f.P.temp << "\n"
           << "Kf:     " << f.P.Kf << "\n"
           << "conc:   " << f.P.conc << "\n"
           << "height: " << f.P.height << "\n";

    myfile << "\n";
    myfile << "Integration parameters used: \n";
    myfile << "\n";
    myfile << "dt:       " << dt << "\n"
           << "T:        " << total_time << "\n"
           << "eps:		   " << tol << "\n"
           << "tSave:    " << tSave << "\n"
           << "no. non-integrated: "
           << f.mask.raw().rows() - f.mask.raw().cast<size_t>().sum() << "\n";
    myfile.close();

  } else
    cout << "Unable to open file";
}

void Integrator::getStatusLog(std::string nameOfFile, std::size_t frame,
                              double areaError, double volumeError,
                              double bendingError, double faceError,
                              std::string inputMesh) {
  ofstream myfile(nameOfFile);
  if (myfile.is_open()) {
    myfile << "Input Mesh: " << inputMesh << "\n";
    myfile << "Final parameter: \n";
    myfile << "\n";
    myfile << "Kb(bare):     " << f.P.Kb << "\n"
           << "Kb(coated):   " << f.P.Kbc << "\n"
           << "H0:     " << f.P.H0 << "\n"
           << "Kse:    " << f.P.Kse << "\n"
           << "Ksl:    " << f.P.Ksl << "\n"
           << "Kst:    " << f.P.Kst << "\n"
           << "Ksg:    " << f.P.Ksg << "\n"
           << "Kv:     " << f.P.Kv << "\n"
           << "gamma:  " << f.P.gamma << "\n"
           << "Vt:     " << f.P.Vt << "\n"
           << "kt:     " << f.P.temp << "\n"
           << "Kf:   " << f.P.Kf << "\n"
           << "conc:   " << f.P.conc << "\n";

    myfile << "\n";
    myfile << "Integration: \n";
    myfile << "\n";
    myfile << "dt:    " << dt << "\n"
           << "T:     " << f.time << "\n"
           << "Frame: " << frame << "\n";

    myfile << "\n";
    myfile << "States: \n";
    myfile << "\n";
    myfile << "Bending Energy:   " << f.E.BE << "\n"
           << "Surface Energy:   " << f.E.sE << "\n"
           << "Pressure Work:    " << f.E.pE << "\n"
           << "Kinetic Work:    " << f.E.kE << "\n"
           << "Chemical Energy:  " << f.E.cE << "\n"
           << "Line tension Energy:  " << f.E.lE << "\n"
           << "Total Energy:     " << f.E.totalE << "\n"
           << "L1 error norm:    " << f.L1ErrorNorm << "\n"
           << "Volume:           " << f.volume << " = "
           << f.volume / f.refVolume << " reduced volume"
           << "\n"
           << "Surface area:     " << f.surfaceArea << " = "
           << f.surfaceArea / f.refSurfaceArea << " target surface area"
           << "\n"
           << "COM (x, y, z):		 "
           << gc::EigenMap<double, 3>(f.vpg->inputVertexPositions)
                      .colwise()
                      .sum() /
                  f.vpg->inputVertexPositions.raw().rows()
           << "\n";

    myfile << "\n";
    myfile << "Errors: \n";
    myfile << "\n";
    myfile << "Bending error:       " << bendingError * 100 << "%"
           << "\n"
           << "Volume error:        " << volumeError * 100 << "%"
           << "\n"
           << "Surface area error:  " << areaError * 100 << "%"
           << "\n"
           << "Face area error:     " << faceError * 100 << "%"
           << "\n";

    myfile << "\n";
    myfile << "Options: \n";
    myfile << "\n";
    myfile << "Is considering protein: " << f.O.isProtein << "\n"
           << "Is vertex shift: " << f.O.isVertexShift << "\n";

    myfile.close();
  } else
    cout << "Unable to open file";
}

} // namespace mem3dg
