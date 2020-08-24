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

#include <iostream>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "mem3dg/solver/ddgsolver.h"
#include "mem3dg/solver/force.h"
#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/typetraits.h"
#include "mem3dg/solver/util.h"

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int viewer(std::string fileName) {
  std::cout << "Initializing the mesh and geometry ...";
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
      gcs::RichSurfaceMeshData::readMeshAndData(fileName);
  ptrVpg = ptrRichData->getGeometry();
  std::cout << "Finished!" << std::endl;

  std::cout << "Reading the vertex properties ...";
  gcs::VertexData<double> meanCurvature =
      ptrRichData->getVertexProperty<double>("mean_curvature");
  gcs::VertexData<double> sponCurvature =
      ptrRichData->getVertexProperty<double>("spon_curvature");
  gcs::VertexData<double> extPressure =
      ptrRichData->getVertexProperty<double>("external_pressure");
  gcs::VertexData<double> physicalPressure =
      ptrRichData->getVertexProperty<double>("physical_pressure");
  gcs::VertexData<double> surfaceTension =
      ptrRichData->getVertexProperty<double>("surface_tension");
  gcs::VertexData<double> bendingPressure =
      ptrRichData->getVertexProperty<double>("bending_pressure");
  /*gcs::VertexData<gc::Vector3> normalForce =
  ptrRichData->getVertexProperty<gc::Vector3>("normal_force");
  gcs::VertexData<gc::Vector3> tangentialForce =
  ptrRichData->getVertexProperty<gc::Vector3>("tangential_force");*/
  std::cout << "Finished!" << std::endl;

  Eigen::Matrix<double, Eigen::Dynamic, 1> meanCurvature_e =
      meanCurvature.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sponCurvature_e =
      sponCurvature.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> extPressure_e = extPressure.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> physicalPressure_e = physicalPressure.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> surfaceTension_e =
      surfaceTension.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> bendingPressure_e = bendingPressure.raw();
  /*Eigen::Matrix<double, Eigen::Dynamic, 3> normalForce_e =
  gc::EigenMap<double, 3>(normalForce); Eigen::Matrix<double,
  Eigen::Dynamic, 3> tangentialForce_e = gc::EigenMap<double,
  3>(tangentialForce);*/

  std::cout << "Opening Polyscope GUI ...";
  polyscope::init();
  polyscope::registerSurfaceMesh("Vesicle surface",
                                 ptrVpg->inputVertexPositions,
                                 ptrMesh->getFaceVertexList());

  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("mean_curvature", meanCurvature_e);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("spon_curvature", sponCurvature_e);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("applied_pressure", extPressure_e);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("surface_tension", surfaceTension_e);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("physical_pressure", physicalPressure_e);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("bending_pressure", bendingPressure_e);
  /*polyscope::getSurfaceMesh("Vesicle
  surface")->addVertexVectorQuantity("tangential_force", tangentialForce_e);
  polyscope::getSurfaceMesh("Vesicle
  surface")->addVertexVectorQuantity("normal_force", normalForce_e);*/
  std::cout << "Finished!" << std::endl;
  polyscope::show();
  
  return 0;
}

int genIcosphere(size_t nSub, std::string path, double R) {
  std::cout << "Constructing " << nSub << "-subdivided icosphere of radius "
            << R << " ...";

  /// initialize mesh and vpg
  std::unique_ptr<gcs::HalfedgeMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;

  /// initialize icosphere
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;
  ddgsolver::icosphere(coords, polygons, nSub, R);
  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);
  // writeSurfaceMesh(*ptrMesh, *ptrVpg, path);
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);
  richData.write(path);

  std::cout << "Finished!" << std::endl;
  return 0;
}

int driver(std::string inputMesh, std::string refMesh, bool isTuftedLaplacian,
           double mollifyFactor, bool isVertexShift, double Kb, double H0,
           double Kse, double Kst, double Ksl, std::vector<double> Ksg, 
           std::vector<double>Kv, double Vt,
           double gamma, double kt, size_t ptInd, double kf, double conc,
           double height, double radius, double h, double T, double eps,
           double closeZone, double increment, double tSave, double tMollify,
           std::string outputDir) {

  /// physical parameters
  double sigma = sqrt(2 * gamma * kt / h);
  ddgsolver::Parameters p{Kb, H0, Ksg[0], Kst, Ksl, Kse,  Kv[0], gamma, Vt,
                          kt, sigma, ptInd, kf,  conc, height, radius};

  std::cout << "Loading input mesh " << inputMesh << " ...";
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  /*std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
  gcs::RichSurfaceMeshData::readMeshAndData(inputMesh); <- this returns no
  connectivity for UVsphere.ply ptrVpg = ptrRichData->getGeometry();*/
  std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(inputMesh);
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);
  std::cout << "Finished!" << std::endl;

  std::cout << "Loading reference mesh " << refMesh << " ...";
  std::unique_ptr<gcs::SurfaceMesh> ptrRefMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrRefVpg = ptrVpg->copy();
  std::cout << "Finished!" << std::endl;

  std::cout << "Initiating the system ...";
  ddgsolver::Force f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p,
                     isTuftedLaplacian, mollifyFactor, isVertexShift);
  std::cout << "Finished!" << std::endl;

  std::cout << "Solving the system ..." << std::endl;
  ddgsolver::integration::velocityVerlet(f, h, T, eps, closeZone, increment,
                                         Kv[1], Ksg[1],
                                         tSave, tMollify, inputMesh, outputDir);

  return 0;
}
