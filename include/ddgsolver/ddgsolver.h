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

#pragma once
#include <geometrycentral/surface/surface_mesh.h>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int viewer(std::string fileName);

int genIcosphere(size_t nSub, std::string path, double R);

int driver(std::string inputMesh, std::string refMesh, bool isTuftedLaplacian, 
           double mollifyFactor, bool isVertexShift, double Kb, double H0,
           double Kse, double Ksl, double Ksg, double Kv, double Vt,
           double gamma, double kt, size_t ptInd, double Kf, double conc, double height, double radius,
           double h, double T, double eps, double closeZone, double increment,
           double tSave, double tMollify, std::string outputDir);
