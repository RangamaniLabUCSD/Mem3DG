
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

#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/integrator.h"
#include "ddgsolver/util.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int main() {
  std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;

  /// initialize icosphere
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  ddgsolver::tetrahedron(coords, polygons);

  // auto makeNormedVertex = [](double x, double y, double z) -> gc::Vector3 {
  //   return gc::Vector3{std::move(x), std::move(y), std::move(z)}.normalize();
  // };

  // coords.emplace_back(makeNormedVertex(1, 0, 0));
  // coords.emplace_back(makeNormedVertex(-1, 0, 0));
  // coords.emplace_back(makeNormedVertex(0, 1, 0));
  
  // polygons.emplace_back(std::vector<std::size_t>{0,1,2});

  gc::PolygonSoupMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  std::tie(ptrmesh, ptrvpg) =
      gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

  ddgsolver::Parameters p;
  p.Kb = 0.01;    // Kb
  p.H0 = 1.5;     // H0
  p.Kse = 0.1;    // Kse
  p.Ksl = 1;      // Ksl
  p.Ksg = 2;      // Ksg
  p.Kv = 10;      // Kv
  p.gamma = 1;    // gamma
  p.Vt = 1 * 0.7; // Vt
  p.kt = 0.00001; // Kt

  auto &mesh = *ptrmesh;
  auto &vpg = *ptrvpg;
  ddgsolver::Force f(mesh, vpg, p);
  velocityVerlet(f, 0.005, 0.01, 1e-9);

  // f.getDampingForces();
  // f.getStochasticForces();
  f.getDPDForces();

  std::cout << "Damping Forces: " << std::endl;
  std::cout << ddgsolver::EigenMap<double, 3>(f.dampingForces) << std::endl;

  std::cout << "Stoch. Forces: " << std::endl;
  std::cout << ddgsolver::EigenMap<double, 3>(f.stochasticForces) << std::endl;
  // polyscope::init();
  // polyscope::registerSurfaceMesh("mymesh", ptrvpg->inputVertexPositions,
  //                                ptrmesh->getFaceVertexList());
  // polyscope::show();

  return 0;
} // namespace gcs=::geometrycentral::surfaceintmain()
