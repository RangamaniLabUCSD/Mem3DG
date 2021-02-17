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
namespace integration {

void checkParameters(std::string integrator, System &f, double h) {
  if (integrator == "velocity verlet") {
    if (abs(f.P.sigma -
         sqrt(2 * f.P.gamma * mem3dg::constants::kBoltzmann * f.P.temp / h)) /
            sqrt(2 * f.P.gamma * mem3dg::constants::kBoltzmann * f.P.temp / h) >
        1e-6) {
      throw std::runtime_error(
          "sigma for DPD is not consistent, probably not initialized!");
    }
  }
  if (integrator == "euler") {
    if (f.P.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for euler integration!");
    }
  }
  if (integrator == "conjugate gradient") {
    if (f.P.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for euler integration!");
    }
  }
}

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

void getParameterLog(System &f, double dt, double finalTime, double tolerance,
                     double tSave, std::string inputMesh,
                     std::string outputDir) {
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
           << "T:        " << finalTime << "\n"
           << "eps:		   " << tolerance << "\n"
           << "tSave:    " << tSave << "\n"
           << "no. non-integrated: "
           << f.mask.raw().rows() - f.mask.raw().cast<size_t>().sum() << "\n";
    myfile.close();

  } else
    cout << "Unable to open file";
}

void getStatusLog(std::string nameOfFile, System &f, double dt, double time,
                  std::size_t frame, double areaError, double volumeError,
                  double bendingError, double faceError, double bendingEnergy,
                  double surfaceEnergy, double pressureEnergy,
                  double kineticEnergy, double chemicalEnergy,
                  double lineEnergy, double totalEnergy, double L2ErrorNorm,
                  bool isTuftedLaplacian, bool isProtein, bool isVertexShift,
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
           << "T:     " << time << "\n"
           << "Frame: " << frame << "\n";

    myfile << "\n";
    myfile << "States: \n";
    myfile << "\n";
    myfile << "Bending Energy:   " << bendingEnergy << "\n"
           << "Surface Energy:   " << surfaceEnergy << "\n"
           << "Pressure Work:    " << pressureEnergy << "\n"
           << "Kinetic Work:    " << kineticEnergy << "\n"
           << "Chemical Energy:  " << chemicalEnergy << "\n"
           << "Line tension Energy:  " << lineEnergy << "\n"
           << "Total Energy:     " << totalEnergy << "\n"
           << "L2 error norm:    " << L2ErrorNorm << "\n"
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
    myfile << "Is tufted laplacian:    " << isTuftedLaplacian << "\n"
           << "Is considering protein: " << isProtein << "\n"
           << "Is vertex shift: " << isVertexShift << "\n";

    myfile.close();
  } else
    cout << "Unable to open file";
}

void getEnergyLog(double time, double bendingEnergy, double surfaceEnergy,
                  double pressureEnergy, double kineticEnergy,
                  double chemicalEnergy, double totalEnergy,
                  std::string outputDir) {
  ofstream myfile(outputDir + "energy.txt", std::ios::app);
  if (myfile.is_open()) {
    myfile << time << "," << bendingEnergy << "," << surfaceEnergy << ","
           << pressureEnergy << "," << kineticEnergy << "," << chemicalEnergy
           << "," << totalEnergy << "\n";
    myfile.close();
  } else
    cout << "Unable to open file";
}

} // namespace integration
} // namespace mem3dg
