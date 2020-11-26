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

int viewPly(std::string fileName);

int viewer(std::string fileName, const bool mean_curvature = 0,
           const bool spon_curvature = 0, const bool ext_pressure = 0,
           const bool physical_pressure = 0, const bool capillary_pressure = 0,
           const bool bending_pressure = 0, const bool line_pressure = 0);

int view_animation(std::string &filename, const bool ref_coord = 0,
                   const bool velocity = 0, const bool mean_curvature = 0,
                   const bool spon_curvature = 0, const bool ext_pressure = 0,
                   const bool physical_pressure = 0,
                   const bool capillary_pressure = 0,
                   const bool bending_pressure = 0,
                   const bool line_pressure = 0, const bool mask = 0,
                   const bool H_H0 = 0);

int genIcosphere(size_t nSub, std::string path, double R);

int driver_ply(const size_t verbosity, std::string inputMesh,
               std::string refMesh, size_t nSub, bool isTuftedLaplacian,
               bool isProtein, double mollifyFactor, bool isVertexShift,
               double Kb, double H0, double sharpness, std::vector<double> r_H0, double Kse,
               double Kst, double Ksl, std::vector<double> Ksg,
               std::vector<double> Kv, double eta, double epsilon, double Bc,
               double Vt, double gamma, double kt, std::vector<double> pt,
               double Kf, double conc, double height, double radius, double h,
               double T, double eps, double closeZone, double increment,
               double tSave, double tMollify, std::string outputDir, double errorJumpLim);

int driver_nc(const size_t verbosity, std::string trajFile,
              std::size_t startingFrame, bool isTuftedLaplacian, bool isProtein,
              double mollifyFactor, bool isVertexShift, double Kb, double H0,
              double sharpness, std::vector<double> r_H0, double Kse, double Kst, double Ksl,
              std::vector<double> Ksg, std::vector<double> Kv, double eta,
              double epsilon, double Bc, double Vt, double gamma, double kt,
              std::vector<double> pt, double Kf, double conc, double height,
              double radius, double h, double T, double eps, double closeZone,
              double increment, double tSave, double tMollify,
              std::string outputDir, double errorJumpLim);
