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
  std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
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
    std::tie(ptrmesh, ptrvpg) =
        gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);
  }
};

TEST_F(ForceCalculationTest, ConsistentForcesTest) {
  ddgsolver::Force f(*ptrmesh, *ptrvpg, p);

  f.getConservativeForces();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> bendingForces1 =
      f.bendingForces.raw();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> pressureForces1 =
      f.pressureForces.raw();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> stretchingForces1 =
      f.stretchingForces.raw();

  f.getConservativeForces();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> bendingForces2 =
      f.bendingForces.raw();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> pressureForces2 =
      f.pressureForces.raw();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> stretchingForces2 =
      f.stretchingForces.raw();

  ASSERT_EQ(bendingForces1, bendingForces2);
  ASSERT_EQ(pressureForces1, pressureForces2);
  ASSERT_EQ(stretchingForces1, stretchingForces2);
};

TEST_F(ForceCalculationTest, OnePassVsReferenceForce) {
  ddgsolver::Force f(*ptrmesh, *ptrvpg, p);

  f.getConservativeForces();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> bendingForces1 =
      f.bendingForces.raw();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> pressureForces1 =
      f.pressureForces.raw();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> stretchingForces1 =
      f.stretchingForces.raw();

  Eigen::Matrix<double, Eigen::Dynamic, 3> bf1 =
    EigenMap<double,3>(f.bendingForces);

  f.getBendingForces();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> bendingForces2 =
      f.bendingForces.raw();

  Eigen::Matrix<double, Eigen::Dynamic, 3> bf2 =
    EigenMap<double, 3>(f.bendingForces);

  f.getPressureForces();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> pressureForces2 =
      f.pressureForces.raw();

  f.getStretchingForces();
  Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> stretchingForces2 =
    f.stretchingForces.raw();

  ASSERT_EQ(bendingForces1, bendingForces2);
  ASSERT_EQ(pressureForces1, pressureForces2);
  ASSERT_EQ(stretchingForces1, stretchingForces2);
};
} // namespace ddgsolver
