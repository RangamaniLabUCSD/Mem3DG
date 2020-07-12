#include <iostream>

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>

#include "ddgsolver/force.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"
#include "ddgsolver/integrator.h"
#include "ddgsolver/icosphere.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int main() {
  std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;

  /// initialize icosphere
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  ddgsolver::icosphere(coords, polygons, 2);
  // ddgsolver::tetrahedron(coords, polygons);

  // auto makeNormedVertex = [](double x, double y, double z) -> gc::Vector3 {
  //   return gc::Vector3{std::move(x), std::move(y), std::move(z)}.normalize();
  // };

  // coords.emplace_back(makeNormedVertex(1, 0, 0));
  // coords.emplace_back(makeNormedVertex(-1, 0, 0));
  // coords.emplace_back(makeNormedVertex(0, 1, 0));
  
  // polygons.emplace_back(std::vector<std::size_t>{0,1,2});

  gcs::PolygonSoupMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  std::tie(ptrmesh, ptrvpg) =
      gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);

  ddgsolver::Parameters p;
  p.Kb = 0;    // Kb
  p.H0 = 0;     // H0
  p.Kse = 0;    // Kse
  p.Ksl = 0;      // Ksl
  p.Ksg = 0;      // Ksg
  p.Kv = 0;      // Kv
  p.gamma = 5;    // gamma
  p.Vt = 1; // Vt
  p.kt = 0.414;// Kt

  auto &mesh = *ptrmesh;
  auto &vpg = *ptrvpg;
  ddgsolver::Force f(mesh, vpg, p);
  ddgsolver::integration::velocityVerlet(f, 0.005, 0.4, 1e-9, 0.1, "output-file/");

  // std::cout << "Damping Forces: " << std::endl;
  // std::cout << ddgsolver::EigenMap<double, 3>(f.dampingForces) << std::endl;

  // std::cout << "Stoch. Forces: " << std::endl;
  // std::cout << ddgsolver::EigenMap<double, 3>(f.stochasticForces) << std::endl;

  // polyscope::init();
  // polyscope::registerSurfaceMesh("mymesh", ptrvpg->inputVertexPositions,
  //                                ptrmesh->getFaceVertexList());
  // polyscope::show();

  return 0;
} // namespace gcs=::geometrycentral::surfaceintmain()
