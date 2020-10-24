#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/util.h"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

// We are writing 2D data, a 6 x 12 grid
constexpr int nx = 6;
constexpr int ny = 12;

// Return this in event of a problem
constexpr int nc_err = 2;

int main() {
  namespace gc = ::geometrycentral;
  namespace gcs = ::geometrycentral::surface;

  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  ddgsolver::icosphere(coords, polygons, 1, 1);

  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();

  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(soup.polygons, soup.vertexCoordinates);
  // std::unique_ptr<gcs::ManifoldSurfaceMesh> mMesh = ptrMesh->toManifoldMesh();
  ddgsolver::subdivide(ptrMesh, ptrVpg, 2);
  std::cout << ptrMesh->nVertices() << std::endl;
  std::cout << ptrVpg->inputVertexPositions.raw().size() << std::endl;
}
