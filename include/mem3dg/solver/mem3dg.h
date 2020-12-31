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

/**
 * @brief Visualize .ply file in polysope with options of additional quantities
 *
 */
int viewer_ply(std::string fileName, const bool mean_curvature = 0,
           const bool spon_curvature = 0, const bool ext_pressure = 0,
           const bool physical_pressure = 0, const bool capillary_pressure = 0,
           const bool bending_pressure = 0, const bool line_pressure = 0);

/**
 * @brief Run single simulation starting with .ply files
 *
 */
int driver_ply(const size_t verbosity, std::string inputMesh,
               std::string refMesh, size_t nSub, bool isTuftedLaplacian,
               bool isProtein, double mollifyFactor, bool isVertexShift,
               double Kb, double H0, double sharpness, std::vector<double> r_H0,
               double Kse, double Kst, double Ksl, double Ksg, double Kv,
               double eta, double epsilon, double Bc, double Vt, double gamma,
               double temp, std::vector<double> pt, double Kf, double conc,
               double height, double radius, double h, double T, double eps,
               double tSave, std::string outputDir,
               std::string integrationMethod, bool isBacktrack, double rho,
               double c1, double ctol, bool isAugmentedLagrangian);

/**
 * @brief Run forward sweep simulation starting with .ply
 * files
 *
 */
int forwardsweep_ply(std::string inputMesh, std::string refMesh, size_t nSub,
                     bool isTuftedLaplacian, bool isProtein,
                     double mollifyFactor, bool isVertexShift, double Kb,
                     std::vector<double> H0, double sharpness,
                     std::vector<double> r_H0, double Kse, double Kst,
                     double Ksl, double Ksg, double Kv, double eta,
                     double epsilon, double Bc, std::vector<double> Vt,
                     double gamma, double temp, std::vector<double> pt,
                     double Kf, double conc, double height, double radius,
                     double h, double T, double eps, double tSave,
                     std::string outputDir, bool isBacktrack, double rho,
                     double c1, double ctol, bool isAugmentedLagrangian);

#ifdef MEM3DG_WITH_NETCDF
/**
 * @brief Visualize netcdf file in single frame
 *
 */
int snapshot_nc(std::string &filename, int frame, bool isShow, bool isSave,
                std::string screenshotName, const bool ref_coord,
                const bool velocity, const bool mean_curvature,
                const bool spon_curvature, const bool ext_pressure,
                const bool physical_pressure, const bool capillary_pressure,
                const bool bending_pressure, const bool line_pressure,
                const bool mask, const bool H_H0);

/**
 * @brief Animate netcdf file with options of additional quantities
 *
 */
int animation_nc(std::string &filename, const bool ref_coord = 0,
                   const bool velocity = 0, const bool mean_curvature = 0,
                   const bool spon_curvature = 0, const bool ext_pressure = 0,
                   const bool physical_pressure = 0,
                   const bool capillary_pressure = 0,
                   const bool bending_pressure = 0,
                   const bool line_pressure = 0, const bool mask = 0,
                   const bool H_H0 = 0);

/**
 * @brief Run single simulation starting with netcdf files
 *
 */
int driver_nc(const size_t verbosity, std::string trajFile,
              std::size_t startingFrame, bool isTuftedLaplacian, bool isProtein,
              double mollifyFactor, bool isVertexShift, double Kb, double H0,
              double sharpness, std::vector<double> r_H0, double Kse,
              double Kst, double Ksl, double Ksg, double Kv, double eta,
              double epsilon, double Bc, double Vt, double gamma, double temp,
              std::vector<double> pt, double Kf, double conc, double height,
              double radius, double h, double T, double eps, double tSave,
              std::string outputDir, std::string integrationMethod,
              bool isBacktrack, double rho, double c1, double ctol,
              bool isAugmentedLagrangian);

/**
 * @brief Run forward sweep simulation starting with netcdf
 * files
 *
 */
int forwardsweep_nc(std::string trajFile, std::size_t startingFrame,
                    bool isTuftedLaplacian, bool isProtein,
                    double mollifyFactor, bool isVertexShift, double Kb,
                    std::vector<double> H0, double sharpness,
                    std::vector<double> r_H0, double Kse, double Kst,
                    double Ksl, double Ksg, double Kv, double eta,
                    double epsilon, double Bc, std::vector<double> Vt,
                    double gamma, double temp, std::vector<double> pt,
                    double Kf, double conc, double height, double radius,
                    double h, double T, double eps, double tSave,
                    std::string outputDir, bool isBacktrack, double rho,
                    double c1, double ctol, bool isAugmentedLagrangian);

#endif