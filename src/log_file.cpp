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

#include "mem3dg/solver/force.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/version.h"

#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <fstream>
#include <iostream>
using namespace std;

namespace ddgsolver {
namespace integration {

void getParameterLog(Force &f, double dt,
                     double finalTime, double tolerance,
                     double tSave, std::string inputMesh,
                     std::string outputDir) {
  ofstream myfile(outputDir + "parameter.txt");
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
           << "kt:     " << f.P.kt << "\n"
           << "sigma:  " << f.P.sigma << "\n"
           << "ptInd:  " << f.P.ptInd << "\n"
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
           << f.mask.rows() - f.mask.cast<size_t>().sum() << "\n";
    myfile.close();

  } else
    cout << "Unable to open file";
}

void getStatusLog(std::string nameOfFile, Force &f, double dt, double time, std::size_t frame, double areaError,
                   double volumeError, double bendingError, double faceError, double bendingEnergy, double surfaceEnergy, 
                   double pressureEnergy, double chemicalEnergy, double totalEnergy,
                  double L2ErrorNorm, bool isTuftedLaplacian, bool isProtein, bool isVertexShift,
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
           << "kt:     " << f.P.kt << "\n"
           << "sigma:  " << f.P.sigma << "\n"
           << "ptInd:  " << f.P.ptInd << "\n"
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
           << "Chemical Energy:  " << chemicalEnergy << "\n"
           << "Total Energy:     " << totalEnergy << "\n" 
           << "L2 error norm:    " << L2ErrorNorm << "\n"
           << "Volume:           " << f.volume << " = "
           << f.volume / f.refVolume << " reduced volume"
           << "\n"
           << "Surface area:     " << f.surfaceArea << " = "
           << f.surfaceArea / f.targetSurfaceArea << " target surface area"
           << "\n"
           << "COM (x, y, z):		 "
           << gc::EigenMap<double, 3>(f.vpg.inputVertexPositions).colwise().sum() /
                  f.vpg.inputVertexPositions.raw().rows()
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

} // namespace integration
} // namespace ddgsolver
