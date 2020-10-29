#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

#include "mem3dg/solver/force.h"
#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/util.h"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

// We are writing 2D data, a 6 x 12 grid
constexpr int nx = 6;
constexpr int ny = 12;

// Return this in event of a problem
constexpr int nc_err = 2;

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

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

int main() {

  std::string inputMesh = "./icosphere.ply";
  std::string refMesh = inputMesh;
  genIcosphere(1, inputMesh, 1);
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::cout << "Loading input mesh " << inputMesh << " ...";
  std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(inputMesh);
  std::cout << "Finished!" << std::endl;

  std::cout << "Loading reference mesh " << refMesh << " ...";
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrRefMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrRefVpg_;
  std::tie(ptrRefMesh, ptrRefVpg_) = gcs::readManifoldSurfaceMesh(refMesh);
  std::cout << "Finished!" << std::endl;

  std::unique_ptr<gcs::VertexPositionGeometry> ptrRefVpg =
      ptrRefVpg_->reinterpretTo(*ptrMesh);
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  double Kb = 8.22e-5, H0 = 0, sharpness = 10, r_H0 = 100, Kst = 0, Ksl = 0,
         Kse = 0, epsilon = 15e-5, Bc = 40, gamma = 0, Vt = 0.7, Kf = 0,
         conc = 25, height = 0, radius = 100, kt = 0, h = 1e-5, Kv = 5e-2,
         Ksg = 0.1, mollifyFactor = 1e-3;
  size_t ptInd = 1;

  bool isProtein = false, isTuftedLaplacian = false, isVertexShift = false;

  double sigma = sqrt(2 * gamma * kt / h);
  if (ptrMesh->hasBoundary() && (Vt != 1.0)) {
    Vt = 1.0;
    std::cout << "Geometry is a patch, so change Vt to 1.0!" << std::endl;
  }

  std::cout << "Initiating the system ...";
  ddgsolver::Parameters p{Kb,    H0,    sharpness, r_H0, Ksg,    Kst,   Ksl,
                          Kse,   Kv,    epsilon,   Bc,   gamma,  Vt,    kt,
                          sigma, ptInd, Kf,        conc, height, radius};
  ddgsolver::Force f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                     isTuftedLaplacian, mollifyFactor, isVertexShift);
  std::cout << "Finished!" << std::endl;

  // std::cout << "Solving the system ..." << std::endl;
  // ddgsolver::integration::velocityVerlet(f, h, T, eps, closeZone, increment,
  //                                        Kv[1], Ksg[1], tSave, tMollify,
  //                                        inputMesh, outputDir);

  return 0;
}