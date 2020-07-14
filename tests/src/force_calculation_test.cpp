#include <iostream>

#include <gtest/gtest.h>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>

#include <Eigen/Core>

#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/util.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class ForceCalculationTest : public testing::Test {
protected:
  /// initialize mesh and vpg
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  Parameters p;

  ForceCalculationTest() {
    /// physical parameters
    p.Kb = 0.02;    // Kb
    p.H0 = 2;       // H0
    p.Kse = 1;      // Kse
    p.Ksl = 3;      // Ksl
    p.Ksg = 1;      // Ksg
    p.Kv = 2;       // Kv
    p.gamma = 1;    // gamma
    p.Vt = 1 * 0.5; // Vt
    p.kt = 0.0001;  // Kt
    p.ptInd = 1;
    p.extF = 2 * 0;
    p.conc = 25;

    std::vector<gc::Vector3> coords;
    std::vector<std::vector<std::size_t>> polygons;

    icosphere(coords, polygons, 2);

    gcs::PolygonSoupMesh soup(polygons, coords);
    soup.mergeIdenticalVertices();
    std::tie(ptrMesh, ptrVpg) =
        gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);
  }
};

TEST_F(ForceCalculationTest, ConsistentForcesTest) {
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  ddgsolver::Force f(*ptrMesh, *ptrVpg, richData, p);

  f.getConservativeForces();
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> bendingForces1
    = ddgsolver::EigenMap<double, 3>(f.bendingForces);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> pressureForces1
    = ddgsolver::EigenMap<double, 3>(f.pressureForces);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> stretchingForces1
    = ddgsolver::EigenMap<double, 3>(f.stretchingForces);

  f.getConservativeForces();
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> bendingForces2
    = ddgsolver::EigenMap<double, 3>(f.bendingForces);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> pressureForces2
    = ddgsolver::EigenMap<double, 3>(f.pressureForces);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> stretchingForces2
    = ddgsolver::EigenMap<double, 3>(f.stretchingForces);

  ASSERT_TRUE((bendingForces1 - bendingForces2).norm() < 1e-12);
  ASSERT_TRUE((stretchingForces1 - stretchingForces2).norm() < 1e-12);
  ASSERT_TRUE((pressureForces1 - pressureForces2).norm() < 1e-12);
};

TEST_F(ForceCalculationTest, OnePassVsReferenceForce) {
  gcs::RichSurfaceMeshData richData(*ptrMesh);
  ddgsolver::Force f(*ptrMesh, *ptrVpg, richData, p);

  f.getConservativeForces();
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> bendingForces1
    = ddgsolver::EigenMap<double, 3>(f.bendingForces);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> pressureForces1
    = ddgsolver::EigenMap<double, 3>(f.pressureForces);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> stretchingForces1
    = ddgsolver::EigenMap<double, 3>(f.stretchingForces);

  f.getBendingForces();
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> bendingForces2
    = ddgsolver::EigenMap<double, 3>(f.bendingForces);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> pressureForces2
    = ddgsolver::EigenMap<double, 3>(f.pressureForces);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> stretchingForces2
    = ddgsolver::EigenMap<double, 3>(f.stretchingForces);

  ASSERT_TRUE((bendingForces1 - bendingForces2).norm() < 1e-12);
  ASSERT_TRUE((stretchingForces1 - stretchingForces2).norm() < 1e-12);
  ASSERT_TRUE((pressureForces1 - pressureForces2).norm() < 1e-12);
};
} // namespace ddgsolver
