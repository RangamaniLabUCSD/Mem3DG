
#include <iostream>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/vector3.h>

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "ddgsolver/ddgsolver.h"
#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/integrator.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int viewer(std::string fileName) {
  /// initialize mesh and vpg
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
      gcs::RichSurfaceMeshData::readMeshAndData(fileName);
  ptrVpg = ptrRichData->getGeometry();

  gcs::VertexData<double> meanCurvature =
      ptrRichData->getVertexProperty<double>("mean_curvature");
  gcs::VertexData<double> extForce =
      ptrRichData->getVertexProperty<double>("external_force");
  gcs::VertexData<double> normalForce =
      ptrRichData->getVertexProperty<double>("normal_force");
  gcs::VertexData<double> tangentialForce =
      ptrRichData->getVertexProperty<double>("tangential_force");

  Eigen::Matrix<double, Eigen::Dynamic, 1> meanCurvature_e =
      meanCurvature.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> extForce_e = extForce.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> normalForce_e = normalForce.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> tangentialForce_e =
      tangentialForce.raw();

  polyscope::init();
  polyscope::registerSurfaceMesh("Vesicle surface",
                                 ptrVpg->inputVertexPositions,
                                 ptrMesh->getFaceVertexList());

  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("mean_curvature", meanCurvature_e);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("applied_force", extForce_e);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("tangential_force", tangentialForce_e);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("normal_force", normalForce_e);
  polyscope::show();
}

int genIcosphere(size_t nSub, std::string path) {
  /// initialize mesh and vpg
  std::unique_ptr<gcs::HalfedgeMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;

  /// initialize icosphere
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;
  ddgsolver::icosphere(coords, polygons, nSub);
  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);
  // writeSurfaceMesh(*ptrMesh, *ptrVpg, path);
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);
  richData.write(path);

  return 0;
}

int driver(std::string inputMesh, double Kb, double H0, double Kse, double Ksl,
           double Ksg, double Kv, double Vt, double gamma, double kt,
           size_t ptInd, double extF, double conc, double h, double T,
           double eps, double tSave, std::string outputDir) {

  /// physical parameters
  double sigma = sqrt(2 * gamma * kt / h);
  ddgsolver::Parameters p{Kb, H0, Ksl,   Ksg,   Kse,  Kv,  gamma,
                          Vt, kt, sigma, ptInd, extF, conc};

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

  /// run the program based on "run"
  ddgsolver::Force f(*ptrMesh, *ptrVpg, richData, p);
  ddgsolver::integration::velocityVerlet(f, h, T, eps, tSave, outputDir);

  return 0;
}
