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

  physicalForce = (f.mask.raw().cast<double>()).array() *
                  (f.M * (f.bendingPressure.raw() + f.capillaryPressure.raw() +
                          f.externalPressure.raw() + f.insidePressure.raw()) +
                   f.lineCapillaryForce.raw())
                      .array();

  DPDForce = f.mask.raw().cast<double>().array() *
             rowwiseDotProduct((EigenMap<double, 3>(f.dampingForce) +
                                EigenMap<double, 3>(f.stochasticForce)),
                               vertexAngleNormal_e)
                 .array();

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
        f.P.lambdaSG += f.P.Ksg * (f.surfaceArea - f.targetSurfaceArea) /
                        f.targetSurfaceArea;
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
        f.P.lambdaSG += f.P.Ksg * (f.surfaceArea - f.targetSurfaceArea) /
                        f.targetSurfaceArea;
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
    fd.writeRefVolume(f.refVolume);
    fd.writeRefSurfArea(f.targetSurfaceArea);
  }
#endif
}

void Integrator::saveData() {
  // save variable to richData and save ply file
  if (verbosity > 3 && !f.O.isGrowMesh) {
    saveRichData();
    char buffer[50];
    sprintf(buffer, "/frame%d", (int)frame);
    f.richData.write(outputDir + buffer + ".ply");
  }

#ifdef MEM3DG_WITH_NETCDF
  // save variable to netcdf traj file
  if (verbosity > 0 && !f.O.isGrowMesh) {
    saveNetcdfData();
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
              << "E_pot: " << f.E.potE << "\n"
              << "|e|L1: " << f.L1ErrorNorm << "\n"
              << "H: [" << f.H.raw().minCoeff() << "," << f.H.raw().maxCoeff()
              << "]"
              << "\n"
              << "K: [" << f.K.raw().minCoeff() << "," << f.K.raw().maxCoeff()
              << "]" << std::endl;
    // << "COM: "
    // << gc::EigenMap<double,
    // 3>(f.vpg->inputVertexPositions).colwise().sum() /
    //         f.vpg->inputVertexPositions.raw().rows()
    // << "\n"
  }
  // break loop if EXIT flag is on
  if (EXIT) {
    if (verbosity > 0 && !f.O.isGrowMesh) {
      std::cout << "Simulation " << (SUCCESS ? "finished" : "failed")
                << ", and data saved to " + outputDir << std::endl;
      if (verbosity > 2) {
        saveRichData();
        f.richData.write(outputDir + "/out.ply");
      }
    }
  }
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

void Integrator::saveRichData() {
  gcs::VertexData<double> fn(*f.mesh);
  fn.fromVector(f.M_inv * physicalForce);
  gcs::VertexData<double> fl(*f.mesh);
  fl.fromVector(f.M_inv * f.lineCapillaryForce.raw());

  f.richData.addVertexProperty("mean_curvature", f.H);
  f.richData.addVertexProperty("gauss_curvature", f.K);
  f.richData.addVertexProperty("spon_curvature", f.H0);
  f.richData.addVertexProperty("external_pressure", f.externalPressure);
  f.richData.addVertexProperty("physical_pressure", fn);
  f.richData.addVertexProperty("capillary_pressure", f.capillaryPressure);
  f.richData.addVertexProperty("bending_pressure", f.bendingPressure);
  f.richData.addVertexProperty("line_tension_pressure", fl);
}

void Integrator::saveRichData(std::string plyName) {
  saveRichData();
  f.richData.write(outputDir + plyName);
}

#ifdef MEM3DG_WITH_NETCDF
void Integrator::saveNetcdfData() {
  frame = fd.getNextFrameIndex();

  // write time
  fd.writeTime(frame, f.time);

  // write geometry
  fd.writeCoords(frame, EigenMap<double, 3>(f.vpg->inputVertexPositions));
  fd.writeTopoFrame(frame, getFaceVertexMatrix(*f.mesh));
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
  if (f.O.isProtein) {
    fd.writeProteinDensity(frame, f.proteinDensity.raw());
  }

  // write pressures
  fd.writeBendingPressure(frame, f.bendingPressure.raw());
  fd.writeCapillaryPressure(frame, f.capillaryPressure.raw());
  fd.writeLinePressure(frame, f.M_inv * f.lineCapillaryForce.raw());
  fd.writeInsidePressure(frame, f.insidePressure.raw());
  fd.writeExternalPressure(frame, f.externalPressure.raw());
  fd.writePhysicalPressure(frame, f.M_inv * physicalForce);

  // write energies
  fd.writeBendEnergy(frame, f.E.BE);
  fd.writeSurfEnergy(frame, f.E.sE);
  fd.writePressEnergy(frame, f.E.pE);
  fd.writeKineEnergy(frame, f.E.kE);
  fd.writeChemEnergy(frame, f.E.cE);
  fd.writeLineEnergy(frame, f.E.lE);
  fd.writeTotalEnergy(frame, f.E.totalE);

  // write Norms
  fd.writeL1ErrorNorm(frame, f.L1ErrorNorm);
  fd.writeL1BendNorm(frame, f.computeL1Norm(f.M * f.bendingPressure.raw()));
  fd.writeL1SurfNorm(frame, f.computeL1Norm(f.M * f.capillaryPressure.raw()));
  fd.writeL1PressNorm(frame, f.computeL1Norm(f.M * f.insidePressure.raw()));
  fd.writeL1LineNorm(frame, f.computeL1Norm(f.lineCapillaryForce.raw()));
}
#endif

void Integrator::getParameterLog(std::string inputMesh) {
  ofstream myfile(outputDir + "/parameter.txt");
  if (myfile.is_open()) {
    myfile << "Mem3DG Version: " << MEM3DG_VERSION << "\n";
    myfile << "Input Mesh:     " << inputMesh << "\n";
    myfile << "Physical parameters used: \n";
    myfile << "\n";
    myfile << "Kb:     " << f.P.Kb << "\n"
           << "H0:     " << f.P.H0 << "\n"
           << "Kse:    " << f.P.Kse << "\n"
           << "Ksl:    " << f.P.Ksl << "\n"
           << "Kst:    " << f.P.Kst << "\n"
           << "Ksg:    " << f.P.Ksg << "\n"
           << "Kv:     " << f.P.Kv << "\n"
           << "gamma:  " << f.P.gamma << "\n"
           << "Vt:     " << f.P.Vt << "\n"
           << "kt:     " << f.P.temp << "\n"
           << "sigma:  " << f.P.sigma << "\n"
           << "ptInd:  " << f.theVertex.getIndex() << "\n"
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
    myfile << "Kb:     " << f.P.Kb << "\n"
           << "H0:     " << f.P.H0 << "\n"
           << "Kse:    " << f.P.Kse << "\n"
           << "Ksl:    " << f.P.Ksl << "\n"
           << "Kst:    " << f.P.Kst << "\n"
           << "Ksg:    " << f.P.Ksg << "\n"
           << "Kv:     " << f.P.Kv << "\n"
           << "gamma:  " << f.P.gamma << "\n"
           << "Vt:     " << f.P.Vt << "\n"
           << "kt:     " << f.P.temp << "\n"
           << "sigma:  " << f.P.sigma << "\n"
           << "ptInd:  " << f.theVertex.getIndex() << "\n"
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
           << f.surfaceArea / f.targetSurfaceArea << " target surface area"
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
