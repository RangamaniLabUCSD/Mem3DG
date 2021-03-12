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
#include "mem3dg/mem3dg"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <geometrycentral/surface/surface_mesh.h>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

/**
 * @brief Run single simulation starting with .ply files
 *
 */

std::unique_ptr<mem3dg::System> system_ply(
    const size_t verbosity, std::string inputMesh, std::string refMesh,
    size_t nSub, bool isReducedVolume, bool isProtein, bool isLocalCurvature,
    bool isVertexShift, bool isEdgeFlip, bool isGrowMesh, double Kb, double H0,
    std::vector<double> r_H0, double Kse, double Kst, double Ksl, double Ksg,
    double Kv, double eta, double epsilon, double Bc, double Vt, double cam,
    double gamma, double temp, std::vector<double> pt, double Kf, double conc,
    double height, double radius, double h, double T, double eps, double tSave,
    std::string outputDir, std::string integrationMethod, bool isBacktrack,
    double rho, double c1, double ctol, bool isAugmentedLagrangian);

/**
 * @brief Run single simulation starting with .ply files
 *
 */
int driver_ply(const size_t verbosity, std::string inputMesh,
               std::string refMesh, size_t nSub, bool isReducedVolume,
               bool isProtein, bool isLocalCurvature, bool isVertexShift,
               bool isEdgeFlip, bool isGrowMesh, double Kb, double H0,
               std::vector<double> r_H0, double Kse, double Kst, double Ksl,
               double Ksg, double Kv, double eta, double epsilon, double Bc,
               double Vt, double cam, double gamma, double temp,
               std::vector<double> pt, double Kf, double conc, double height,
               double radius, double h, double T, double eps, double tSave,
               std::string outputDir, std::string integrationMethod,
               bool isBacktrack, double rho, double c1, double ctol,
               bool isAugmentedLagrangian, bool isAdaptiveStep);

/**
 * @brief Run forward sweep simulation starting with .ply
 * files
 *
 */
int forwardsweep_ply(std::string inputMesh, std::string refMesh, size_t nSub,
                     bool isReducedVolume, bool isProtein,
                     bool isLocalCurvature, bool isVertexShift, bool isEdgeFlip,
                     bool isGrowMesh, double Kb, std::vector<double> H0,
                     std::vector<double> r_H0, double Kse, double Kst,
                     double Ksl, double Ksg, double Kv, double eta,
                     double epsilon, double Bc, std::vector<double> Vt,
                     std::vector<double> cam, double gamma, double temp,
                     std::vector<double> pt, double Kf, double conc,
                     double height, double radius, double h, double T,
                     double eps, double tSave, std::string outputDir,
                     bool isBacktrack, double rho, double c1, double ctol,
                     bool isAugmentedLagrangian, bool isAdaptiveStep);

#ifdef MEM3DG_WITH_NETCDF
/**
 * @brief Run single simulation starting with netcdf files
 *
 */
int driver_nc(const size_t verbosity, std::string trajFile, int startingFrame,
              int nSub, bool isContinue, bool isReducedVolume, bool isProtein,
              bool isLocalCurvature, bool isVertexShift, bool isEdgeFlip,
              bool isGrowMesh, double Kb, double H0, std::vector<double> r_H0,
              double Kse, double Kst, double Ksl, double Ksg, double Kv,
              double eta, double epsilon, double Bc, double Vt, double cam,
              double gamma, double temp, std::vector<double> pt, double Kf,
              double conc, double height, double radius, double h, double T,
              double eps, double tSave, std::string outputDir,
              std::string integrationMethod, bool isBacktrack, double rho,
              double c1, double ctol, bool isAugmentedLagrangian,
              bool isAdaptiveStep);

/**
 * @brief Run forward sweep simulation starting with netcdf
 * files
 *
 */
int forwardsweep_nc(std::string trajFile, int startingFrame, int nSub,
                    bool isContinue, bool isReducedVolume, bool isProtein,
                    bool isLocalCurvature, bool isVertexShift, bool isEdgeFlip,
                    bool isGrowMesh, double Kb, std::vector<double> H0,
                    std::vector<double> r_H0, double Kse, double Kst,
                    double Ksl, double Ksg, double Kv, double eta,
                    double epsilon, double Bc, std::vector<double> Vt,
                    std::vector<double> cam, double gamma, double temp,
                    std::vector<double> pt, double Kf, double conc,
                    double height, double radius, double h, double T,
                    double eps, double tSave, std::string outputDir,
                    bool isBacktrack, double rho, double c1, double ctol,
                    bool isAugmentedLagrangian, bool isAdaptiveStep);

#endif