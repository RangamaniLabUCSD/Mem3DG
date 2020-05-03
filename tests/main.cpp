
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

  icosphere(coords, polygons, 2);

  gc::PolygonSoupMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();

  std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
  std::tie(ptrmesh, ptrvpg) =
      gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

  auto &mesh = *ptrmesh;
  auto &vpg = *ptrvpg;
  double h = 0.01;
  // initiate force object f
  Force f(mesh, vpg, h);

  f.getBendingForces(1.0, 0);
  f.getStretchingForces(1.0, 1.0);
  f.getPressureForces(1.0, 0.7);
  f.getDampingForces(1.0);
  f.getStochasticForces(1.0);

  for (gcs::Vertex v : mesh.vertices()) {
    std::cout << "bending force" << f.bendingForces[v] << std::endl;
    std::cout << "stretching force" << f.stretchingForces[v] << std::endl;
    std::cout << "pressure force" << f.pressureForces[v] << std::endl;
    std::cout << "damping force" << f.dampingForces[v] << std::endl;
    std::cout << "stochastic force" << f.stochasticForces[v] << std::endl;
  }

  polyscope::init();
  polyscope::registerSurfaceMesh("myMesh", ptrvpg->inputVertexPositions,
                                 ptrmesh->getFaceVertexList());
  polyscope::show();

  return 0;
}
