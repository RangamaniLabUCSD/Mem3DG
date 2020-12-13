#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

#include "mem3dg/solver/ddgsolver.h"
#include "mem3dg/solver/force.h"
#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/util.h"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int main() {
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;
  ddgsolver::icosphere(coords, polygons, 1, 1);
  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(soup.polygons,
      soup.vertexCoordinates);
  gcs::VertexPositionGeometry *ptrRefVpg =
      new gcs::VertexPositionGeometry(*ptrMesh,
      ptrVpg->inputVertexPositions);

  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  double Kb = 8.22e-5, H0 = 0, sharpness = 10, Kst = 10, Ksl = 0,
         Kse = 0, epsilon = 15e-5, Bc = 40, gamma = 0, Vt = 0.7, Kf = 0,
         conc = 25, height = 0, radius = 0.9, kt = 0, h = 5e-4, Kv = 5e-2, eta = 0,
         Ksg = 0.1, mollifyFactor = 1e-3;
  std::vector<double> pt = {1,1,1};
  std::vector<double> r_H0 = {100, 100};

  bool isProtein = false, isTuftedLaplacian = false, isVertexShift = false;

  double sigma = sqrt(2 * gamma * kt / h);
  if (ptrMesh->hasBoundary() && (Vt != 1.0)) {
    Vt = 1.0;
    std::cout << "Geometry is a patch, so change Vt to 1.0!" << std::endl;
  }

  std::cout << "Initiating the system ...";
  ddgsolver::Parameters p{Kb,    H0,    sharpness, r_H0, Ksg,    Kst,   Ksl,
                          Kse,   Kv, eta,   epsilon,   Bc,   gamma,  Vt,    kt,
                          sigma + 1e-18, pt, Kf, conc, height, radius};
  ddgsolver::Force f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                     isTuftedLaplacian, mollifyFactor, isVertexShift);
  std::cout << "Finished!" << std::endl;

  std::cout << "Solving the system ..." << std::endl;
  double T = 3, eps = 0.002, closeZone = 1000, increment = 0, tSave = 1e-1,
  tMollify = 100; size_t verbosity = 0;
  // ddgsolver::integration::velocityVerlet(f, h, T, eps, closeZone, increment,
  //                                        Kv, Ksg, tSave, tMollify,
  //                                        verbosity);
  ddgsolver::integration::conjugateGradient(
        f, h, T, eps, closeZone, increment, Kv, Ksg, tSave, tMollify,
        verbosity);
  delete ptrRefVpg;
  return 0;
}