

#include <iostream>

#include <gtest/gtest.h>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "ddgsolver/icosphere.h"
#include "ddgsolver/util.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class MappingUtilityTest : public testing::Test {};

TEST_F(MappingUtilityTest, MappingTest) {
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

  auto pos1 = mapVecToEigen<double, 3>(vpg.inputVertexPositions);
  auto pos2 = mapVecToEigen<double, 3>(vpg.inputVertexPositions.rawdata());
  auto pos3 = mapVecToEigen<gc::Vector3, 1>(vpg.inputVertexPositions);
  auto pos4 = mapVecToEigen<gc::Vector3, 1>(vpg.inputVertexPositions.rawdata());

  for (std::size_t i = 0; i < mesh.nVertices(); ++i) {
    std::size_t xidx = 3 * i;
    std::size_t yidx = 3 * i + 1;
    std::size_t zidx = 3 * i + 2;
    // check x values
    ASSERT_EQ(vpg.inputVertexPositions[i][0], pos1(xidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][0], pos2(xidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][0], pos3(i)[0]);
    ASSERT_EQ(vpg.inputVertexPositions[i][0], pos4(i)[0]);

    // check y values
    ASSERT_EQ(vpg.inputVertexPositions[i][1], pos1(yidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][1], pos2(yidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][1], pos3(i)[1]);
    ASSERT_EQ(vpg.inputVertexPositions[i][1], pos4(i)[1]);

    // check z values
    ASSERT_EQ(vpg.inputVertexPositions[i][2], pos1(zidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][2], pos2(zidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][2], pos3(i)[2]);
    ASSERT_EQ(vpg.inputVertexPositions[i][2], pos4(i)[2]);

    vpg.inputVertexPositions[i][0] += 1;
    ASSERT_EQ(vpg.inputVertexPositions[i][0], pos1(xidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][0], pos2(xidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][0], pos3(i)[0]);
    ASSERT_EQ(vpg.inputVertexPositions[i][0], pos4(i)[0]);

    pos1(yidx) -= 3.14;
    ASSERT_EQ(vpg.inputVertexPositions[i][1], pos1(yidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][1], pos2(yidx));
    ASSERT_EQ(vpg.inputVertexPositions[i][1], pos3(i)[1]);
    ASSERT_EQ(vpg.inputVertexPositions[i][1], pos4(i)[1]);
  }
}