
#include <iostream>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

// overload << to print vector;
template <typename T>
std::ostream &operator<<(std::ostream &output, const std::vector<T> &v) {
  output << "[";
  for (size_t i = 0; i != v.size() - 1; ++i) {
    output << v[i] << ",";
  }
  output << v[v.size() - 1];
  output << "]";
  return output;
}

int main() {

  // initialization of the geometry and mesh
  // Construct an icosahedron
  // Construct a std::vector of Vector 3
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  icosphere(coords, polygons, 5);

  gc::PolygonSoupMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();

  // std::cout << soup.polygons << std::endl;

  std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
  std::tie(ptrmesh, ptrvpg) =
      gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

  auto &mesh = *ptrmesh;
  auto &vpg = *ptrvpg;
  double h = 0.01;

  // // initiate force object f
  // Force f(mesh, vpg, h);

  // // f.pcg_test();

  // std::cout << "Sizeof Force: " << sizeof(f) << std::endl;

  // // calculate the bending force
  // f.bendingForces(1.0, 0);
  // std::cout << "bending force" << f.bendingForces << std::endl;

  // // calculate the stretching force
  // f.stretchingForces(1.0, 1.0);
  // std::cout << "stretching force" << f.stretchingForces << std::endl;

  // // calculate the stretching force
  // f.pressure_force(1.0, 0.7);
  // std::cout << "pressure force" << f.pressureForces << std::endl;

  // f.damping_force(1.0);
  // std::cout << "damping force" << f.dampingForces << std::endl;

  // f.stochastic_force(1.0);
  // std::cout << "stochastic force" << f.stochasticForces << std::endl;

  polyscope::init();
  polyscope::registerSurfaceMesh("myMesh", ptrvpg->inputVertexPositions,
                                 ptrmesh->getFaceVertexList());
  polyscope::show();

  return 0;
}
