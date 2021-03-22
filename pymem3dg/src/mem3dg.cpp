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

int driver_ply(const size_t verbosity, std::string inputMesh,
               std::string refMesh, size_t nSub, bool isReducedVolume,
               bool isProtein, bool isLocalCurvature, bool isVertexShift,
               bool isEdgeFlip, bool isGrowMesh, bool isRefMesh,
               bool isFloatVertex, bool isLaplacianMeanCurvature, double Kb,
               double Kbc, double H0, std::vector<double> r_H0, double Kse,
               double Kst, double Ksl, double Ksg, double Kv, double eta,
               double epsilon, double Bc, double Vt, double cam, double gamma,
               double temp, std::vector<double> pt, double Kf, double conc,
               double height, double radius, double h, double T, double eps,
               double tSave, std::string outputDir,
               std::string integrationMethod, bool isBacktrack, double rho,
               double c1, double ctol, bool isAugmentedLagrangian,
               size_t restartNum, bool isAdaptiveStep) {
  /*std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
  gcs::RichSurfaceMeshData::readMeshAndData(inputMesh); <- this returns no
  connectivity for UVsphere.ply ptrVpg = ptrRichData->getGeometry();*/

  /// Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);
  // pybind11::scoped_interpreter guard{};

  /// Initialize parameter struct
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);

  mem3dg::Parameters p{Kb,    Kbc, H0,      r_H0, Ksg,    Kst,   Ksl, Kse,
                       Kv,    eta, epsilon, Bc,   gamma,  Vt,    cam, temp,
                       sigma, pt,  Kf,      conc, height, radius};

  mem3dg::Options o{isVertexShift,    isProtein,     isReducedVolume,
                    isLocalCurvature, isEdgeFlip,    isGrowMesh,
                    isRefMesh,        isFloatVertex, isLaplacianMeanCurvature};

  /// Initialize the system
  mem3dg::System f(inputMesh, refMesh, nSub, p, o);

  /// Time integration / optimization
  if (integrationMethod == "velocity verlet") {
    mem3dg::VelocityVerlet integrator(f, h, isAdaptiveStep, T, tSave, eps,
                                      outputDir, "/traj.nc", verbosity);
    integrator.integrate();
  } else if (integrationMethod == "euler") {
    mem3dg::Euler integrator(f, h, isAdaptiveStep, T, tSave, eps, outputDir,
                             "/traj.nc", verbosity, isBacktrack, rho, c1);
    integrator.integrate();
  } else if (integrationMethod == "conjugate gradient") {
    mem3dg::ConjugateGradient integrator(
        f, h, isAdaptiveStep, T, tSave, eps, outputDir, "/traj.nc", verbosity,
        isBacktrack, rho, c1, ctol, isAugmentedLagrangian, restartNum);
    integrator.integrate();
  } else if (integrationMethod == "BFGS") {
    mem3dg::BFGS integrator(f, h, isAdaptiveStep, T, tSave, eps, outputDir,
                            "/traj.nc", verbosity, isBacktrack, rho, c1, ctol,
                            isAugmentedLagrangian);
    integrator.integrate();
  }

  return 0;
}

int forwardsweep_ply(
    std::string inputMesh, std::string refMesh, size_t nSub,
    bool isReducedVolume, bool isProtein, bool isLocalCurvature,
    bool isVertexShift, bool isEdgeFlip, bool isGrowMesh, bool isRefMesh,
    bool isFloatVertex, bool isLaplacianMeanCurvature, double Kb, double Kbc,
    std::vector<double> H0, std::vector<double> r_H0, double Kse, double Kst,
    double Ksl, double Ksg, double Kv, double eta, double epsilon, double Bc,
    std::vector<double> Vt, std::vector<double> cam, double gamma, double temp,
    std::vector<double> pt, double Kf, double conc, double height,
    double radius, double h, double T, double eps, double tSave,
    std::string outputDir, bool isBacktrack, double rho, double c1, double ctol,
    bool isAugmentedLagrangian, size_t restartNum, bool isAdaptiveStep) {

  /// Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);

  /// Initialize parameter struct
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  mem3dg::Parameters p{Kb,    Kbc, H0[0],   r_H0, Ksg,    Kst,   Ksl,    Kse,
                       Kv,    eta, epsilon, Bc,   gamma,  Vt[0], cam[0], temp,
                       sigma, pt,  Kf,      conc, height, radius};

  mem3dg::Options o{isVertexShift,    isProtein,     isReducedVolume,
                    isLocalCurvature, isEdgeFlip,    isGrowMesh,
                    isRefMesh,        isFloatVertex, isLaplacianMeanCurvature};

  /// Initialize the system
  mem3dg::System f(inputMesh, refMesh, nSub, p, o);

  /// Time integration / optimization
  mem3dg::FeedForwardSweep integrator(
      f, h, isAdaptiveStep, T, tSave, eps, outputDir, "/traj.nc", 3,
      isBacktrack, rho, c1, ctol, isAugmentedLagrangian, restartNum, H0,
      (isReducedVolume) ? Vt : cam);
  integrator.sweep();

  return 0;
}

#ifdef MEM3DG_WITH_NETCDF
int driver_nc(const size_t verbosity, std::string trajFile, int startingFrame,
              int nSub, bool isContinue, bool isReducedVolume, bool isProtein,
              bool isLocalCurvature, bool isVertexShift, bool isEdgeFlip,
              bool isGrowMesh, bool isRefMesh, bool isFloatVertex,
              bool isLaplacianMeanCurvature, double Kb, double Kbc, double H0,
              std::vector<double> r_H0, double Kse, double Kst, double Ksl,
              double Ksg, double Kv, double eta, double epsilon, double Bc,
              double Vt, double cam, double gamma, double temp,
              std::vector<double> pt, double Kf, double conc, double height,
              double radius, double h, double T, double eps, double tSave,
              std::string outputDir, std::string integrationMethod,
              bool isBacktrack, double rho, double c1, double ctol,
              bool isAugmentedLagrangian, size_t restartNum,
              bool isAdaptiveStep) {

  // Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);
  // pybind11::scoped_interpreter guard{};

  // Initialize parameter struct
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  mem3dg::Parameters p{Kb,    Kbc, H0,      r_H0, Ksg,    Kst,   Ksl, Kse,
                       Kv,    eta, epsilon, Bc,   gamma,  Vt,    cam, temp,
                       sigma, pt,  Kf,      conc, height, radius};

  mem3dg::Options o{isVertexShift,    isProtein,     isReducedVolume,
                    isLocalCurvature, isEdgeFlip,    isGrowMesh,
                    isRefMesh,        isFloatVertex, isLaplacianMeanCurvature};

  // Initialize the system
  mem3dg::System f(trajFile, startingFrame, nSub, isContinue, p, o);

  /// Time integration / optimization
  if (integrationMethod == "velocity verlet") {
    mem3dg::VelocityVerlet integrator(f, h, isAdaptiveStep, T, tSave, eps,
                                      outputDir, "/traj.nc", verbosity);
    integrator.integrate();
  } else if (integrationMethod == "euler") {
    mem3dg::Euler integrator(f, h, isAdaptiveStep, T, tSave, eps, outputDir,
                             "/traj.nc", verbosity, isBacktrack, rho, c1);
    integrator.integrate();
  } else if (integrationMethod == "conjugate gradient") {
    mem3dg::ConjugateGradient integrator(
        f, h, isAdaptiveStep, T, tSave, eps, outputDir, "/traj.nc", verbosity,
        isBacktrack, rho, c1, ctol, isAugmentedLagrangian, restartNum);
    integrator.integrate();
  } else if (integrationMethod == "BFGS") {
    mem3dg::BFGS integrator(f, h, isAdaptiveStep, T, tSave, eps, outputDir,
                            "/traj.nc", verbosity, isBacktrack, rho, c1, ctol,
                            isAugmentedLagrangian);
    integrator.integrate();
  }

  return 0;
}

int forwardsweep_nc(
    std::string trajFile, int startingFrame, int nSub, bool isContinue,
    bool isReducedVolume, bool isProtein, bool isLocalCurvature,
    bool isVertexShift, bool isEdgeFlip, bool isGrowMesh, bool isRefMesh,
    bool isFloatVertex, bool isLaplacianMeanCurvature, double Kb, double Kbc,
    std::vector<double> H0, std::vector<double> r_H0, double Kse, double Kst,
    double Ksl, double Ksg, double Kv, double eta, double epsilon, double Bc,
    std::vector<double> Vt, std::vector<double> cam, double gamma, double temp,
    std::vector<double> pt, double Kf, double conc, double height,
    double radius, double h, double T, double eps, double tSave,
    std::string outputDir, bool isBacktrack, double rho, double c1, double ctol,
    bool isAugmentedLagrangian, size_t restartNum, bool isAdaptiveStep) {

  /// Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);
  // pybind11::scoped_interpreter guard{};

  /// Initialize parameter struct
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  mem3dg::Parameters p{Kb,    Kbc, H0[0],   r_H0, Ksg,    Kst,   Ksl,    Kse,
                       Kv,    eta, epsilon, Bc,   gamma,  Vt[0], cam[0], temp,
                       sigma, pt,  Kf,      conc, height, radius};

  mem3dg::Options o{isVertexShift,    isProtein,     isReducedVolume,
                    isLocalCurvature, isEdgeFlip,    isGrowMesh,
                    isRefMesh,        isFloatVertex, isLaplacianMeanCurvature};

  /// Initialize the system
  mem3dg::System f(trajFile, startingFrame, nSub, isContinue, p, o);

  /// Time integration / optimization
  mem3dg::FeedForwardSweep integrator(
      f, h, isAdaptiveStep, T, tSave, eps, outputDir, "/traj.nc", 3,
      isBacktrack, rho, c1, ctol, isAugmentedLagrangian, restartNum, H0,
      (isReducedVolume) ? Vt : cam);
  integrator.sweep();

  return 0;
}

#endif
