
#include <Eigen/Core>
#include <iostream>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/ply_halfedge_mesh_data.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "ddgsolver/icosphere.h"
#include "ddgsolver/util.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int main() {

  Eigen::Matrix<double, 10, 3> foo;

  for(int i = 0; i < 10; ++i){
    foo(i,0) = 3*i;
    foo(i,1) = 3*i+1;
    foo(i,2) = 3*i+2;
  }

  std::cout << foo << std::endl;

  Eigen::Matrix<double, 1, 3> bar;
  bar << 0, 3, 6;

  std::cout << bar << std::endl;

  foo = foo.rowwise() - bar;

  std::cout << foo << std::endl; 



  // std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
  // std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;

  // /// initialize icosphere
  // std::vector<gc::Vector3> coords;
  // std::vector<std::vector<std::size_t>> polygons;
  // ddgsolver::tetrahedron(coords, polygons);
  // gc::PolygonSoupMesh soup(polygons, coords);
  // soup.mergeIdenticalVertices();
  // std::tie(ptrmesh, ptrvpg) =
  //     gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);



  // polyscope::init();
  // polyscope::registerSurfaceMesh("mymesh", ptrvpg->inputVertexPositions,
  //                                ptrmesh->getFaceVertexList());
  // polyscope::show();

  return 0;
} // namespace gcs=::geometrycentral::surfaceintmain()
