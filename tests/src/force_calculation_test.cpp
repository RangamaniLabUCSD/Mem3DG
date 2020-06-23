#include <iostream>

#include <gtest/gtest.h>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/ply_halfedge_mesh_data.h>
#include <geometrycentral/utilities/vector3.h>

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
  ddgsolver::Parameters p;

  ForceCalculationTest() {
    /// physical parameters 
    p.Kb = 0.02;			//Kb
    p.H0 = 2;				//H0
    p.Kse = 0;      //Kse
    p.Ksl = 3;				//Ksl
    p.Ksg = 0;				//Ksg
    p.Kv = 1;			  //Kv
    p.gamma = 1;				//gamma
    p.Vt = 1 * 0.73;			//Vt
    p.kt = 0.0001;		//Kt 
    p.ptInd = 1;
    p.extF = 2 * 0;
    p.conc = 25;

    std::vector<gc::Vector3> coords;
    std::vector<std::vector<std::size_t>> polygons;

    icosphere(coords, polygons, 2);

    gc::PolygonSoupMesh soup(polygons, coords);
    soup.mergeIdenticalVertices();
    std::tie(ptrmesh, ptrvpg) = gcs::makeHalfedgeAndGeometry(
      soup.polygons, soup.vertexCoordinates, true);
  }
};

TEST_F(ForceCalculationTest, ConservativeForceTest) {
  auto& mesh = *ptrmesh;
  auto& vpg = *ptrvpg;
  ddgsolver::Force f(mesh, vpg, p);
  f.getConservativeForces();
  auto conservativeForce1 = ddgsolver::EigenMap<double, 3>(f.bendingForces) +
    ddgsolver::EigenMap<double, 3>(f.stretchingForces) +
    ddgsolver::EigenMap<double, 3>(f.pressureForces);
  f.getConservativeForces();
  auto conservativeForce2 = ddgsolver::EigenMap<double, 3>(f.bendingForces) +
    ddgsolver::EigenMap<double, 3>(f.stretchingForces) +
    ddgsolver::EigenMap<double, 3>(f.pressureForces);
  ASSERT_EQ(conservativeForce1, conservativeForce2);
};
}