
#include <iostream>

#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>

#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/integrator.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

double driver(std::string option, int nSub, double H0, double Vt, double h,
              double T, double eps) {
  /// physical parameters
  ddgsolver::Parameters p;
  p.Kb = 0.03;    // Kb
  p.H0 = H0;      // H0
  p.Kse = 0.1;    // Kse
  p.Ksl = 3;      // Ksl
  p.Ksg = 0;      // Ksg
  p.Kv = 2;       // Kv
  p.gamma = 1;    // gamma
  p.Vt = 1 * Vt;  // Vt
  p.kt = 0.00001; // Kt

  //// choose the run
  // std::string run = "visualization"; // 1. "integration 2. "visualization

  ///// Choose the starting mesh
  // std::string option = "sphere"; // 1. "sphere" 2. "continue" 3.
  // "nameOfTheFile" = "output-file/Vt_%d_H0_%d.ply"

  /// initialize mesh and vpg
  std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;

  if (option == "continue") {
    char buffer[50];
    sprintf(buffer, "output-file/Vt_%d_H0_%d.ply", int(p.Vt * 100),
            int(p.H0 * 100));
    std::tie(ptrmesh, ptrvpg) = gcs::loadMesh(buffer);
  } else if (option == "sphere") {
    /// initialize icosphere
    std::vector<gc::Vector3> coords;
    std::vector<std::vector<std::size_t>> polygons;
    ddgsolver::icosphere(coords, polygons, nSub);
    gcs::PolygonSoupMesh soup(polygons, coords);
    soup.mergeIdenticalVertices();
    std::tie(ptrmesh, ptrvpg) = gcs::makeHalfedgeAndGeometry(
        soup.polygons, soup.vertexCoordinates);
  } else {
    std::tie(ptrmesh, ptrvpg) = gcs::loadMesh(option);
  }

  auto &mesh = *ptrmesh;
  auto &vpg = *ptrvpg;

  /// solve
  ddgsolver::Force f(mesh, vpg, p);
  velocityVerlet(f, h, T, eps);

  /// save the .ply file
  gcs::RichSurfaceMeshData data(mesh);
  data.addGeometry(vpg);
  char buffer[50];
  sprintf(buffer, "output-file/Vt_%d_H0_%d.ply", int(p.Vt * 100),
          int(p.H0 * 100));
  data.write(buffer);

  return 0;
}
