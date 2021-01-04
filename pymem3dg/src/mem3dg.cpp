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

void signalHandler(int signum) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";
  exit(signum);
}

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
               double c1, double ctol, bool isAugmentedLagrangian) {
  /*std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
  gcs::RichSurfaceMeshData::readMeshAndData(inputMesh); <- this returns no
  connectivity for UVsphere.ply ptrVpg = ptrRichData->getGeometry();*/

  /// Activate signal handling
  signal(SIGINT, signalHandler);
  // pybind11::scoped_interpreter guard{};

  /// Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrRefMesh_;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrRefVpg_;
  gcs::VertexPositionGeometry *ptrRefVpg;

  /// Load input mesh and geometry
  std::cout << "Loading input mesh " << inputMesh << " ...";
  std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(inputMesh);
  std::cout << "Finished!" << std::endl;

  /// Load input reference mesh and geometry
  std::cout << "Loading reference mesh " << refMesh << " ...";
  std::tie(ptrRefMesh_, ptrRefVpg_) = gcs::readManifoldSurfaceMesh(refMesh);
  std::cout << "Finished!" << std::endl;

  /// Subdivide the mesh and geometry objects
  if (nSub > 0) {
    std::cout << "Subdivide input and reference mesh " << nSub
              << " time(s) ...";
    // mem3dg::subdivide(ptrMesh, ptrVpg, nSub);
    // mem3dg::subdivide(ptrRefMesh, ptrRefVpg, nSub);
    mem3dg::loopSubdivide(ptrMesh, ptrVpg, nSub);
    mem3dg::loopSubdivide(ptrRefMesh_, ptrRefVpg_, nSub);
    std::cout << "Finished!" << std::endl;
  }

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  mem3dg::loadRefMesh(
      ptrMesh, ptrRefVpg,
      gc::EigenMap<double, 3>(ptrRefVpg_->inputVertexPositions));

  /// Initializa richData for ply file
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  /// Initialize parameter struct
  std::cout << "Initializing the system ...";
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  if (ptrMesh->hasBoundary() && Vt != 1.0) {
    throw std::runtime_error("Vt has to be 1 for open boundary simulation!");
  }
  mem3dg::Parameters p{Kb,   H0,    sharpness, r_H0,    Ksg,  Kst,    Ksl,
                       Kse,  Kv,    eta,       epsilon, Bc,   gamma,  Vt,
                       temp, sigma, pt,        Kf,      conc, height, radius};

  /// Initialize the system
  mem3dg::System f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                   isVertexShift, isTuftedLaplacian, mollifyFactor);
  std::cout << "Finished!" << std::endl;

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (integrationMethod == "velocity verlet") {
    mem3dg::integration::velocityVerlet(f, h, 0, T, tSave, eps, verbosity,
                                        outputDir);
  } else if (integrationMethod == "euler") {
    if (p.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for euler integration!");
    }
    mem3dg::integration::euler(f, h, 0, T, tSave, eps, verbosity, outputDir,
                               isBacktrack, rho, c1);
  } else if (integrationMethod == "conjugate gradient") {
    if (p.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for CG optimization!");
    }
    mem3dg::integration::conjugateGradient(
        f, h, 0, T, tSave, eps, ctol, verbosity, outputDir, isBacktrack, rho,
        c1, isAugmentedLagrangian, "/traj.nc");
  }

  /// Delete non unique pointer
  delete ptrRefVpg;

  return 0;
}

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
                     double c1, double ctol, bool isAugmentedLagrangian) {

  /// Activate signal handling
  signal(SIGINT, signalHandler);

  /// Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrRefMesh_;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrRefVpg_;
  gcs::VertexPositionGeometry *ptrRefVpg;

  /// Load input mesh and geometry
  std::cout << "Loading input and reference mesh " << inputMesh << " ...";
  std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(inputMesh);
  std::tie(ptrRefMesh_, ptrRefVpg_) = gcs::readManifoldSurfaceMesh(refMesh);
  std::cout << "Finished!" << std::endl;

  /// Subdivide the mesh and geometry objects
  if (nSub > 0) {
    std::cout << "Subdivide input and reference mesh " << nSub
              << " time(s) ...";
    mem3dg::loopSubdivide(ptrMesh, ptrVpg, nSub);
    mem3dg::loopSubdivide(ptrRefMesh_, ptrRefVpg_, nSub);
    std::cout << "Finished!" << std::endl;
  }

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  mem3dg::loadRefMesh(
      ptrMesh, ptrRefVpg,
      gc::EigenMap<double, 3>(ptrRefVpg_->inputVertexPositions));

  /// Initializa richData for ply file
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  /// Initialize parameter struct
  std::cout << "Initializing the system ...";
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  if (ptrMesh->hasBoundary() && Vt[0] != 1.0) {
    throw std::runtime_error("Vt has to be 1 for open boundary simulation!");
  }
  mem3dg::Parameters p{Kb,   H0[0], sharpness, r_H0,    Ksg,  Kst,    Ksl,
                       Kse,  Kv,    eta,       epsilon, Bc,   gamma,  Vt[0],
                       temp, sigma, pt,        Kf,      conc, height, radius};

  /// Initialize the system
  mem3dg::System f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                   isVertexShift, isTuftedLaplacian, mollifyFactor);
  std::cout << "Finished!" << std::endl;

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (p.gamma != 0) {
    throw std::runtime_error("gamma has to be 0 for CG optimization!");
  }
  mem3dg::integration::feedForwardSweep(f, H0, Vt, h, T, tSave, eps, ctol,
                                        outputDir, isBacktrack, rho, c1,
                                        isAugmentedLagrangian);

  /// Delete non unique pointer
  delete ptrRefVpg;

  return 0;
}

#ifdef MEM3DG_WITH_NETCDF

/**
 * @brief netcdf file frame reader
 *
 * @param fd reference to netcdf trajectory file
 * @param frame reference to the frame index
 *
 */
void getNcFrame(mem3dg::TrajFile &fd, int &frame) {
  int maxFrame = fd.getNextFrameIndex() - 1;
  if (frame > maxFrame || frame < -(maxFrame + 1)) {
    throw std::runtime_error("Snapshot frame exceed limiting frame index!");
  } else if (frame < 0) {
    frame = frame + maxFrame + 1;
  }
}

int driver_nc(const size_t verbosity, std::string trajFile,
              int startingFrame, bool isTuftedLaplacian, bool isProtein,
              double mollifyFactor, bool isVertexShift, double Kb, double H0,
              double sharpness, std::vector<double> r_H0, double Kse,
              double Kst, double Ksl, double Ksg, double Kv, double eta,
              double epsilon, double Bc, double Vt, double gamma, double temp,
              std::vector<double> pt, double Kf, double conc, double height,
              double radius, double h, double T, double eps, double tSave,
              std::string outputDir, std::string integrationMethod,
              bool isBacktrack, double rho, double c1, double ctol,
              bool isAugmentedLagrangian) {

  /// Activate signal handling
  signal(SIGINT, signalHandler);
  // pybind11::scoped_interpreter guard{};

  /// alias eigen matrix
  using EigenVectorX3D =
      Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

  /// Declare variables
  double time;
  EigenVectorX3D coords;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  gcs::VertexPositionGeometry *ptrRefVpg;

  /// Load input mesh and geometry
  std::cout << "Loading input mesh from trajectory file " << trajFile;
  mem3dg::TrajFile fd = mem3dg::TrajFile::openReadOnly(trajFile);
  getNcFrame(fd, startingFrame);
  std::tie(time, coords) = fd.getTimeAndCoords(startingFrame);
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(coords, fd.getTopology());
  std::cout << " and continuing from time t = " << time << " ...";
  std::cout << "Finished!" << std::endl;

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  std::cout << "Loading reference mesh from trajectory file" << trajFile
            << " ...";
  mem3dg::loadRefMesh(ptrMesh, ptrRefVpg, fd.getRefcoordinate());
  std::cout << "Finished!" << std::endl;

  /// Initializa richData for ply file
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  /// Initialize parameter struct
  std::cout << "Initializing the system ...";
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  if (ptrMesh->hasBoundary() && Vt != 1.0) {
    throw std::runtime_error("Vt has to be 1 for open boundary simulation!");
  }
  mem3dg::Parameters p{Kb,   H0,    sharpness, r_H0,    Ksg,  Kst,    Ksl,
                       Kse,  Kv,    eta,       epsilon, Bc,   gamma,  Vt,
                       temp, sigma, pt,        Kf,      conc, height, radius};

  /// Initialize the system
  mem3dg::System f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                   isVertexShift, isTuftedLaplacian, mollifyFactor);
  gc::EigenMap<double, 3>(f.vel) = fd.getVelocity(startingFrame);
  f.proteinDensity.raw() = fd.getProteinDensity(startingFrame);
  f.update_Vertex_positions();
  std::cout << "Finished!" << std::endl;

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (integrationMethod == "velocity verlet") {
    mem3dg::integration::velocityVerlet(f, h, time, T, tSave, eps, verbosity,
                                        outputDir);
  } else if (integrationMethod == "euler") {
    if (p.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for euler integration!");
    }
    mem3dg::integration::euler(f, h, time, T, tSave, eps, verbosity, outputDir,
                               isBacktrack, rho, c1);
  } else if (integrationMethod == "conjugate gradient") {
    if (p.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for CG optimization!");
    }
    mem3dg::integration::conjugateGradient(
        f, h, time, T, tSave, eps, ctol, verbosity, outputDir, isBacktrack, rho,
        c1, isAugmentedLagrangian, "/traj.nc");
  }

  /// Delete non unique pointer
  delete ptrRefVpg;

  return 0;
}

int forwardsweep_nc(std::string trajFile, int startingFrame,
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
                    double c1, double ctol, bool isAugmentedLagrangian) {

  /// Activate signal handling
  signal(SIGINT, signalHandler);
  // pybind11::scoped_interpreter guard{};

  /// alias eigen matrix
  using EigenVectorX3D =
      Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

  /// Declare variables
  double time;
  EigenVectorX3D coords;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  gcs::VertexPositionGeometry *ptrRefVpg;

  /// Load input mesh and geometry
  std::cout << "Loading input mesh from trajectory file " << trajFile;
  mem3dg::TrajFile fd = mem3dg::TrajFile::openReadOnly(trajFile);
  getNcFrame(fd, startingFrame);
  std::tie(time, coords) = fd.getTimeAndCoords(startingFrame);
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(coords, fd.getTopology());
  std::cout << " and continuing from time t = " << time << " ...";
  std::cout << "Finished!" << std::endl;

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  std::cout << "Loading reference mesh from trajectory file" << trajFile
            << " ...";
  mem3dg::loadRefMesh(ptrMesh, ptrRefVpg, fd.getRefcoordinate());
  std::cout << "Finished!" << std::endl;

  /// Initializa richData for ply file
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  /// Initialize parameter struct
  std::cout << "Initializing the system ...";
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);
  if (ptrMesh->hasBoundary() && Vt[0] != 1.0) {
    throw std::runtime_error("Vt has to be 1 for open boundary simulation!");
  }
  mem3dg::Parameters p{Kb,   H0[0], sharpness, r_H0,    Ksg,  Kst,    Ksl,
                       Kse,  Kv,    eta,       epsilon, Bc,   gamma,  Vt[0],
                       temp, sigma, pt,        Kf,      conc, height, radius};

  /// Initialize the system
  mem3dg::System f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                   isVertexShift, isTuftedLaplacian, mollifyFactor);
  gc::EigenMap<double, 3>(f.vel) = fd.getVelocity(startingFrame);
  std::cout << "Finished!" << std::endl;

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (p.gamma != 0) {
    throw std::runtime_error("gamma has to be 0 for CG optimization!");
  }
  mem3dg::integration::feedForwardSweep(f, H0, Vt, h, T, tSave, eps, ctol,
                                        outputDir, isBacktrack, rho, c1,
                                        isAugmentedLagrangian);

  /// Delete non unique pointer
  delete ptrRefVpg;

  return 0;
}

#endif
