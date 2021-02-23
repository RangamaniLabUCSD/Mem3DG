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
#include "mem3dg/solver/system.h"
#include "mem3dg/solver/version.h"

#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <fstream>
#include <iostream>
using namespace std;

namespace mem3dg {

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
  if (verbosity > 3) {
    saveRichData();
    char buffer[50];
    sprintf(buffer, "/frame%d", (int)frame);
    f.richData.write(outputDir + buffer + ".ply");
  }

#ifdef MEM3DG_WITH_NETCDF
  // save variable to netcdf traj file
  if (verbosity > 0) {
    saveNetcdfData(frame, fd);
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
    if (verbosity > 0) {
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
void Integrator::saveNetcdfData(size_t &frame, TrajFile &fd) {

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
  fd.writeTime(frame, f.time);

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
  fd.writeL1ErrorNorm(frame, f.L1ErrorNorm);
  fd.writeL1BendNorm(frame,
                     f.computeL1Norm(f.M * EigenMap<double, 3>(f.bendingPressure)));
  fd.writeL1SurfNorm(frame,
                     f.computeL1Norm(f.M * EigenMap<double, 3>(f.capillaryPressure)));
  fd.writeL1PressNorm(
      frame, f.computeL1Norm(f.M * f.insidePressure *
                             gc::EigenMap<double, 3>(f.vpg->vertexNormals)));
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
    myfile << "Is considering protein: " << f.isProtein << "\n"
           << "Is vertex shift: " << f.isVertexShift << "\n";

    myfile.close();
  } else
    cout << "Unable to open file";
}

} // namespace mem3dg
