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

#include <csignal>
#include <iostream>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/mem3dg.h"
#include "mem3dg/solver/mesh.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/typetraits.h"
#include "mem3dg/solver/util.h"

#include <pybind11/embed.h>

#include "igl/loop.h"

#include <Eigen/Core>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

std::unique_ptr<mem3dg::System> system_ply(
    const size_t verbosity, std::string inputMesh, std::string refMesh,
    size_t nSub, bool isReducedVolume, bool isProtein, bool isLocalCurvature,
    bool isVertexShift, double Kb, double H0, double sharpness,
    std::vector<double> r_H0, double Kse, double Kst, double Ksl, double Ksg,
    double Kv, double eta, double epsilon, double Bc, double Vt, double cam,
    double gamma, double temp, std::vector<double> pt, double Kf, double conc,
    double height, double radius, double h, double T, double eps, double tSave,
    std::string outputDir, std::string integrationMethod, bool isBacktrack,
    double rho, double c1, double ctol, bool isAugmentedLagrangian) {
  /*std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
  gcs::RichSurfaceMeshData::readMeshAndData(inputMesh); <- this returns no
  connectivity for UVsphere.ply ptrVpg = ptrRichData->getGeometry();*/

  /// Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);
  // pybind11::scoped_interpreter guard{};

  /// Initialize parameter struct
  std::cout << "Initializing the system ..." << std::endl;
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);

  mem3dg::Parameters p{Kb,    H0,  sharpness, r_H0, Ksg,    Kst,   Ksl, Kse,
                       Kv,    eta, epsilon,   Bc,   gamma,  Vt,    cam, temp,
                       sigma, pt,  Kf,        conc, height, radius};

  std::unique_ptr<mem3dg::System> f(
      new mem3dg::System(inputMesh, refMesh, nSub, p, isReducedVolume,
                         isProtein, isLocalCurvature, isVertexShift));
  return f;
}

int driver_ply(const size_t verbosity, std::string inputMesh,
               std::string refMesh, size_t nSub, bool isReducedVolume,
               bool isProtein, bool isLocalCurvature, bool isVertexShift,
               double Kb, double H0, double sharpness, std::vector<double> r_H0,
               double Kse, double Kst, double Ksl, double Ksg, double Kv,
               double eta, double epsilon, double Bc, double Vt, double cam,
               double gamma, double temp, std::vector<double> pt, double Kf,
               double conc, double height, double radius, double h, double T,
               double eps, double tSave, std::string outputDir,
               std::string integrationMethod, bool isBacktrack, double rho,
               double c1, double ctol, bool isAugmentedLagrangian,
               bool isAdaptiveStep) {
  /*std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
  gcs::RichSurfaceMeshData::readMeshAndData(inputMesh); <- this returns no
  connectivity for UVsphere.ply ptrVpg = ptrRichData->getGeometry();*/

  /// Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);
  // pybind11::scoped_interpreter guard{};

  /// Initialize parameter struct
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);

  mem3dg::Parameters p{Kb,    H0,  sharpness, r_H0, Ksg,    Kst,   Ksl, Kse,
                       Kv,    eta, epsilon,   Bc,   gamma,  Vt,    cam, temp,
                       sigma, pt,  Kf,        conc, height, radius};

  /// Initialize the system
  mem3dg::System f(inputMesh, refMesh, nSub, p, isReducedVolume, isProtein,
                   isLocalCurvature, isVertexShift);

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (integrationMethod == "velocity verlet") {
    mem3dg::integration::velocityVerlet(f, h, T, tSave, eps, verbosity,
                                        isAdaptiveStep, outputDir);
  } else if (integrationMethod == "euler") {
    bool success =
        mem3dg::integration::euler(f, h, T, tSave, eps, verbosity, outputDir,
                                   isBacktrack, rho, c1, isAdaptiveStep);
    if (!success) {
      mem3dg::markFileName(outputDir, "/traj.nc", "_failed");
    }
  } else if (integrationMethod == "conjugate gradient") {
    bool success = mem3dg::integration::conjugateGradient(
        f, h, T, tSave, eps, ctol, verbosity, outputDir, isBacktrack, rho, c1,
        isAugmentedLagrangian, isAdaptiveStep, "/traj.nc");
    // mark "failed" is CG returns false
    if (!success) {
      mem3dg::markFileName(outputDir, "/traj.nc", "_failed");
    }
  }

  return 0;
}

int forwardsweep_ply(std::string inputMesh, std::string refMesh, size_t nSub,
                     bool isReducedVolume, bool isProtein,
                     bool isLocalCurvature, bool isVertexShift, double Kb,
                     std::vector<double> H0, double sharpness,
                     std::vector<double> r_H0, double Kse, double Kst,
                     double Ksl, double Ksg, double Kv, double eta,
                     double epsilon, double Bc, std::vector<double> Vt,
                     std::vector<double> cam, double gamma, double temp,
                     std::vector<double> pt, double Kf, double conc,
                     double height, double radius, double h, double T,
                     double eps, double tSave, std::string outputDir,
                     bool isBacktrack, double rho, double c1, double ctol,
                     bool isAugmentedLagrangian, bool isAdaptiveStep) {

  /// Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);

  /// Initialize parameter struct
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  mem3dg::Parameters p{Kb,    H0[0], sharpness, r_H0,  Ksg,     Kst,
                       Ksl,   Kse,   Kv,        eta,   epsilon, Bc,
                       gamma, Vt[0], cam[0],    temp,  sigma,   pt,
                       Kf,    conc,  height,    radius};

  /// Initialize the system
  mem3dg::System f(inputMesh, refMesh, nSub, p, isReducedVolume, isProtein,
                   isLocalCurvature, isVertexShift);

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  mem3dg::integration::feedForwardSweep(
      f, H0, (isReducedVolume) ? Vt : cam, h, T, tSave, eps, ctol, outputDir,
      isBacktrack, rho, c1, isAugmentedLagrangian, isAdaptiveStep);

  return 0;
}

#ifdef MEM3DG_WITH_NETCDF
int driver_nc(const size_t verbosity, std::string trajFile, int startingFrame,
              int nSub, bool isContinue, bool isReducedVolume, bool isProtein,
              bool isLocalCurvature, bool isVertexShift, double Kb, double H0,
              double sharpness, std::vector<double> r_H0, double Kse,
              double Kst, double Ksl, double Ksg, double Kv, double eta,
              double epsilon, double Bc, double Vt, double cam, double gamma,
              double temp, std::vector<double> pt, double Kf, double conc,
              double height, double radius, double h, double T, double eps,
              double tSave, std::string outputDir,
              std::string integrationMethod, bool isBacktrack, double rho,
              double c1, double ctol, bool isAugmentedLagrangian,
              bool isAdaptiveStep) {

  // Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);
  // pybind11::scoped_interpreter guard{};

  // Initialize parameter struct
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  mem3dg::Parameters p{Kb,    H0,  sharpness, r_H0, Ksg,    Kst,   Ksl, Kse,
                       Kv,    eta, epsilon,   Bc,   gamma,  Vt,    cam, temp,
                       sigma, pt,  Kf,        conc, height, radius};

  // Initialize the system
  mem3dg::System f(trajFile, startingFrame, nSub, isContinue, p,
                   isReducedVolume, isProtein, isLocalCurvature, isVertexShift);

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (integrationMethod == "velocity verlet") {
    mem3dg::integration::velocityVerlet(f, h, T, tSave, eps, verbosity,
                                        isAdaptiveStep, outputDir);
  } else if (integrationMethod == "euler") {
    bool success =
        mem3dg::integration::euler(f, h, T, tSave, eps, verbosity, outputDir,
                                   isBacktrack, rho, c1, isAdaptiveStep);
    if (!success) {
      mem3dg::markFileName(outputDir, "/traj.nc", "_failed");
    }
  } else if (integrationMethod == "conjugate gradient") {
    bool success = mem3dg::integration::conjugateGradient(
        f, h, T, tSave, eps, ctol, verbosity, outputDir, isBacktrack, rho, c1,
        isAugmentedLagrangian, isAdaptiveStep, "/traj.nc");

    // mark "failed" is CG returns false
    if (!success) {
      mem3dg::markFileName(outputDir, "/traj.nc", "_failed");
    }
  }

  return 0;
}

int forwardsweep_nc(std::string trajFile, int startingFrame, int nSub,
                    bool isContinue, bool isReducedVolume, bool isProtein,
                    bool isLocalCurvature, bool isVertexShift, double Kb,
                    std::vector<double> H0, double sharpness,
                    std::vector<double> r_H0, double Kse, double Kst,
                    double Ksl, double Ksg, double Kv, double eta,
                    double epsilon, double Bc, std::vector<double> Vt,
                    std::vector<double> cam, double gamma, double temp,
                    std::vector<double> pt, double Kf, double conc,
                    double height, double radius, double h, double T,
                    double eps, double tSave, std::string outputDir,
                    bool isBacktrack, double rho, double c1, double ctol,
                    bool isAugmentedLagrangian, bool isAdaptiveStep) {

  /// Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);
  // pybind11::scoped_interpreter guard{};

  /// Initialize parameter struct
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  mem3dg::Parameters p{Kb,    H0[0], sharpness, r_H0,  Ksg,     Kst,
                       Ksl,   Kse,   Kv,        eta,   epsilon, Bc,
                       gamma, Vt[0], cam[0],    temp,  sigma,   pt,
                       Kf,    conc,  height,    radius};

  /// Initialize the system
  mem3dg::System f(trajFile, startingFrame, nSub, isContinue, p,
                   isReducedVolume, isProtein, isLocalCurvature, isVertexShift);

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  mem3dg::integration::feedForwardSweep(
      f, H0, (isReducedVolume) ? Vt : cam, h, T, tSave, eps, ctol, outputDir,
      isBacktrack, rho, c1, isAugmentedLagrangian, isAdaptiveStep);

  return 0;
}

#endif
