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

#include "mem3dg/solver/ddgsolver.h"
#include "mem3dg/solver/force.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/mesh.h"
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

int viewPly(std::string fileName) {

  signal(SIGINT, signalHandler);

  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(fileName);

  polyscope::init();
  polyscope::registerSurfaceMesh("Vesicle surface",
                                 ptrVpg->inputVertexPositions,
                                 ptrMesh->getFaceVertexList());
  polyscope::show();

  return 0;
}

int viewer(std::string fileName, const bool mean_curvature,
           const bool spon_curvature, const bool ext_pressure,
           const bool physical_pressure, const bool capillary_pressure,
           const bool bending_pressure, const bool line_pressure) {

  signal(SIGINT, signalHandler);

  /// alias eigen matrix
  using EigenVectorX1D = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using EigenVectorX3D =
      Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
  using EigenTopVec =
      Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
      gcs::RichSurfaceMeshData::readMeshAndData(fileName);
  ptrVpg = ptrRichData->getGeometry();

  polyscope::init();
  polyscope::registerSurfaceMesh("Vesicle surface",
                                 ptrVpg->inputVertexPositions,
                                 ptrMesh->getFaceVertexList());

  if (mean_curvature) {
    gcs::VertexData<double> meanCurvature =
        ptrRichData->getVertexProperty<double>("mean_curvature");
    EigenVectorX1D meanCurvature_e = meanCurvature.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("mean_curvature", meanCurvature_e);
  }
  if (spon_curvature) {
    gcs::VertexData<double> sponCurvature =
        ptrRichData->getVertexProperty<double>("spon_curvature");
    EigenVectorX1D sponCurvature_e = sponCurvature.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("spon_curvature", sponCurvature_e);
  }
  if (ext_pressure) {
    gcs::VertexData<double> extPressure =
        ptrRichData->getVertexProperty<double>("external_pressure");
    EigenVectorX1D extPressure_e = extPressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("applied_pressure", extPressure_e);
  }
  if (physical_pressure) {
    gcs::VertexData<double> physicalPressure =
        ptrRichData->getVertexProperty<double>("physical_pressure");
    EigenVectorX1D physicalPressure_e = physicalPressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("physical_pressure", physicalPressure_e);
  }
  if (capillary_pressure) {
    gcs::VertexData<double> capillaryPressure =
        ptrRichData->getVertexProperty<double>("capillary_pressure");
    EigenVectorX1D capillaryPressure_e = capillaryPressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("surface_tension", capillaryPressure_e);
  }
  if (bending_pressure) {
    gcs::VertexData<double> bendingPressure =
        ptrRichData->getVertexProperty<double>("bending_pressure");
    EigenVectorX1D bendingPressure_e = bendingPressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("bending_pressure", bendingPressure_e);
  }
  if (line_pressure) {
    gcs::VertexData<double> linePressure =
        ptrRichData->getVertexProperty<double>("line_tension_pressure");
    EigenVectorX1D linePressure_e = linePressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("line_tension_pressure", linePressure_e);
  }

  /*gcs::VertexData<gc::Vector3> vertexVelocity =
      ptrRichData->getVertexProperty<gc::Vector3>("vertex_velocity");*/
  /*gcs::VertexData<gc::Vector3> normalForce =
  ptrRichData->getVertexProperty<gc::Vector3>("normal_force");
  gcs::VertexData<gc::Vector3> tangentialForce =
  ptrRichData->getVertexProperty<gc::Vector3>("tangential_force");*/
  // EigenVectorX3D vertexVelocity_e =
  //    ddgsolver::EigenMap<double, 3>(vertexVelocity);
  /*EigenVectorX3D normalForce_e =
  gc::EigenMap<double, 3>(normalForce); Eigen::Matrix<double,
  Eigen::Dynamic, 3> tangentialForce_e = gc::EigenMap<double,
  3>(tangentialForce);*/
  /*polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexVectorQuantity("vertexVelocity", vertexVelocity_e);*/
  /*polyscope::getSurfaceMesh("Vesicle
  surface")->addVertexVectorQuantity("tangential_force", tangentialForce_e);
  polyscope::getSurfaceMesh("Vesicle
  surface")->addVertexVectorQuantity("normal_force", normalForce_e);*/

  std::cout << "Finished!" << std::endl;
  polyscope::show();

  return 0;
}

int driver_ply(const size_t verbosity, std::string inputMesh,
               std::string refMesh, size_t nSub, bool isTuftedLaplacian,
               bool isProtein, double mollifyFactor, bool isVertexShift,
               double Kb, double H0, double sharpness, std::vector<double> r_H0,
               double Kse, double Kst, double Ksl, double Ksg, double Kv,
               double eta, double epsilon, double Bc, double Vt, double gamma,
               double kt, std::vector<double> pt, double Kf, double conc,
               double height, double radius, double h, double T, double eps,
               double tSave, std::string outputDir,
               std::string integrationMethod, bool isBacktrack, double rho,
               double c1, bool isAugmentedLagrangian) {
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
    // ddgsolver::subdivide(ptrMesh, ptrVpg, nSub);
    // ddgsolver::subdivide(ptrRefMesh, ptrRefVpg, nSub);
    ddgsolver::loopSubdivide(ptrMesh, ptrVpg, nSub);
    ddgsolver::loopSubdivide(ptrRefMesh_, ptrRefVpg_, nSub);
    std::cout << "Finished!" << std::endl;
  }

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  ddgsolver::loadRefMesh(
      ptrMesh, ptrRefVpg,
      gc::EigenMap<double, 3>(ptrRefVpg_->inputVertexPositions));

  /// Initializa richData for ply file
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  /// Initialize parameter struct
  std::cout << "Initializing the system ...";
  double sigma = sqrt(2 * gamma * kt / h);
  if (ptrMesh->hasBoundary() && Vt != 1.0) {
    throw std::runtime_error("Vt has to be 1 for open boundary simulation!");
  }
  ddgsolver::Parameters p{Kb,  H0,    sharpness, r_H0,    Ksg,  Kst,    Ksl,
                          Kse, Kv,    eta,       epsilon, Bc,   gamma,  Vt,
                          kt,  sigma, pt,        Kf,      conc, height, radius};

  /// Initialize the system
  ddgsolver::System f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                      isTuftedLaplacian, mollifyFactor, isVertexShift);
  std::cout << "Finished!" << std::endl;

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (integrationMethod == "velocity verlet") {
    ddgsolver::integration::velocityVerlet(f, h, 0, T, tSave, eps, verbosity,
                                           outputDir);
  } else if (integrationMethod == "euler") {
    if (p.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for euler integration!");
    }
    ddgsolver::integration::euler(f, h, 0, T, tSave, eps, verbosity, outputDir,
                                  isBacktrack, rho, c1);
  } else if (integrationMethod == "conjugate gradient") {
    if (p.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for CG optimization!");
    }
    ddgsolver::integration::conjugateGradient(
        f, h, 0, T, tSave, eps, 0.005, verbosity, outputDir, isBacktrack, rho,
        c1, isAugmentedLagrangian, "/traj.nc");
  }

  /// Delete non unique pointer
  delete ptrRefVpg;

  return 0;
}

int driver_ply_sweep(std::string inputMesh, std::string refMesh, size_t nSub,
                     bool isTuftedLaplacian, bool isProtein,
                     double mollifyFactor, bool isVertexShift, double Kb,
                     std::vector<double> H0, double sharpness,
                     std::vector<double> r_H0, double Kse, double Kst,
                     double Ksl, double Ksg, double Kv, double eta,
                     double epsilon, double Bc, std::vector<double> Vt,
                     double gamma, double kt, std::vector<double> pt, double Kf,
                     double conc, double height, double radius, double h,
                     double T, double eps, double tSave, std::string outputDir,
                     bool isBacktrack, double rho, double c1,
                     bool isAugmentedLagrangian) {

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
    ddgsolver::loopSubdivide(ptrMesh, ptrVpg, nSub);
    ddgsolver::loopSubdivide(ptrRefMesh_, ptrRefVpg_, nSub);
    std::cout << "Finished!" << std::endl;
  }

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  ddgsolver::loadRefMesh(
      ptrMesh, ptrRefVpg,
      gc::EigenMap<double, 3>(ptrRefVpg_->inputVertexPositions));

  /// Initializa richData for ply file
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  /// Initialize parameter struct
  std::cout << "Initializing the system ...";
  double sigma = sqrt(2 * gamma * kt / h);
  if (ptrMesh->hasBoundary() && Vt[0] != 1.0) {
    throw std::runtime_error("Vt has to be 1 for open boundary simulation!");
  }
  ddgsolver::Parameters p{Kb,  H0[0], sharpness, r_H0,    Ksg,  Kst,    Ksl,
                          Kse, Kv,    eta,       epsilon, Bc,   gamma,  Vt[0],
                          kt,  sigma, pt,        Kf,      conc, height, radius};

  /// Initialize the system
  ddgsolver::System f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                      isTuftedLaplacian, mollifyFactor, isVertexShift);
  std::cout << "Finished!" << std::endl;

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (p.gamma != 0) {
    throw std::runtime_error("gamma has to be 0 for CG optimization!");
  }
  ddgsolver::integration::feedForwardSweep(f, H0, Vt, h, T, tSave, eps, 0.005,
                                           outputDir, isBacktrack, rho, c1,
                                           isAugmentedLagrangian);

  /// Delete non unique pointer
  delete ptrRefVpg;

  return 0;
}

#ifdef MEM3DG_WITH_NETCDF
int driver_nc(const size_t verbosity, std::string trajFile,
              std::size_t startingFrame, bool isTuftedLaplacian, bool isProtein,
              double mollifyFactor, bool isVertexShift, double Kb, double H0,
              double sharpness, std::vector<double> r_H0, double Kse,
              double Kst, double Ksl, double Ksg, double Kv, double eta,
              double epsilon, double Bc, double Vt, double gamma, double kt,
              std::vector<double> pt, double Kf, double conc, double height,
              double radius, double h, double T, double eps, double tSave,
              std::string outputDir, std::string integrationMethod,
              bool isBacktrack, double rho, double c1,
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
  ddgsolver::TrajFile fd = ddgsolver::TrajFile::openReadOnly(trajFile);
  std::tie(time, coords) = fd.getTimeAndCoords(startingFrame);
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(coords, fd.getTopology());
  std::cout << " and continuing from time t = " << time << " ...";
  std::cout << "Finished!" << std::endl;

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  std::cout << "Loading reference mesh from trajectory file" << trajFile
            << " ...";
  ddgsolver::loadRefMesh(ptrMesh, ptrRefVpg, fd.getRefcoordinate());
  std::cout << "Finished!" << std::endl;

  /// Initializa richData for ply file
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  /// Initialize parameter struct
  std::cout << "Initializing the system ...";
  double sigma = sqrt(2 * gamma * kt / h);
  if (ptrMesh->hasBoundary() && Vt != 1.0) {
    throw std::runtime_error("Vt has to be 1 for open boundary simulation!");
  }
  ddgsolver::Parameters p{Kb,  H0,    sharpness, r_H0,    Ksg,  Kst,    Ksl,
                          Kse, Kv,    eta,       epsilon, Bc,   gamma,  Vt,
                          kt,  sigma, pt,        Kf,      conc, height, radius};

  /// Initialize the system
  ddgsolver::System f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                      isTuftedLaplacian, mollifyFactor, isVertexShift);
  gc::EigenMap<double, 3>(f.vel) = fd.getVelocity(startingFrame);
  std::cout << "Finished!" << std::endl;

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (integrationMethod == "velocity verlet") {
    ddgsolver::integration::velocityVerlet(f, h, time, T, tSave, eps, verbosity,
                                           outputDir);
  } else if (integrationMethod == "euler") {
    if (p.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for euler integration!");
    }
    ddgsolver::integration::euler(f, h, time, T, tSave, eps, verbosity,
                                  outputDir, isBacktrack, rho, c1);
  } else if (integrationMethod == "conjugate gradient") {
    if (p.gamma != 0) {
      throw std::runtime_error("gamma has to be 0 for CG optimization!");
    }
    ddgsolver::integration::conjugateGradient(
        f, h, time, T, tSave, eps, 0.005, verbosity, outputDir, isBacktrack,
        rho, c1, isAugmentedLagrangian, "/traj.nc");
  }

  /// Delete non unique pointer
  delete ptrRefVpg;

  return 0;
}

int driver_nc_sweep(std::string trajFile, std::size_t startingFrame,
                    bool isTuftedLaplacian, bool isProtein,
                    double mollifyFactor, bool isVertexShift, double Kb,
                    std::vector<double> H0, double sharpness,
                    std::vector<double> r_H0, double Kse, double Kst,
                    double Ksl, double Ksg, double Kv, double eta,
                    double epsilon, double Bc, std::vector<double> Vt,
                    double gamma, double kt, std::vector<double> pt, double Kf,
                    double conc, double height, double radius, double h,
                    double T, double eps, double tSave, std::string outputDir,
                    bool isBacktrack, double rho, double c1,
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
  ddgsolver::TrajFile fd = ddgsolver::TrajFile::openReadOnly(trajFile);
  std::tie(time, coords) = fd.getTimeAndCoords(startingFrame);
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(coords, fd.getTopology());
  std::cout << " and continuing from time t = " << time << " ...";
  std::cout << "Finished!" << std::endl;

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  std::cout << "Loading reference mesh from trajectory file" << trajFile
            << " ...";
  ddgsolver::loadRefMesh(ptrMesh, ptrRefVpg, fd.getRefcoordinate());
  std::cout << "Finished!" << std::endl;

  /// Initializa richData for ply file
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  /// Initialize parameter struct
  std::cout << "Initializing the system ...";
  double sigma = sqrt(2 * gamma * kt / h);
  if (ptrMesh->hasBoundary() && Vt[0] != 1.0) {
    throw std::runtime_error("Vt has to be 1 for open boundary simulation!");
  }
  ddgsolver::Parameters p{Kb,  H0[0], sharpness, r_H0,    Ksg,  Kst,    Ksl,
                          Kse, Kv,    eta,       epsilon, Bc,   gamma,  Vt[0],
                          kt,  sigma, pt,        Kf,      conc, height, radius};

  /// Initialize the system
  ddgsolver::System f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                      isTuftedLaplacian, mollifyFactor, isVertexShift);
  gc::EigenMap<double, 3>(f.vel) = fd.getVelocity(startingFrame);
  std::cout << "Finished!" << std::endl;

  /// Time integration / optimization
  std::cout << "Solving the system and saving to " << outputDir << std::endl;
  if (p.gamma != 0) {
    throw std::runtime_error("gamma has to be 0 for CG optimization!");
  }
  ddgsolver::integration::feedForwardSweep(f, H0, Vt, h, T, tSave, eps, 0.005,
                                           outputDir, isBacktrack, rho, c1,
                                           isAugmentedLagrangian);

  /// Delete non unique pointer
  delete ptrRefVpg;

  return 0;
}

#endif
