

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

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class MappingUtilityTest : public testing::Test {
protected:
  MappingUtilityTest() {
    std::vector<gc::Vector3> coords;
    std::vector<std::vector<std::size_t>> polygons;

    icosphere(coords, polygons, 2);

    gc::PolygonSoupMesh soup(polygons, coords);
    soup.mergeIdenticalVertices();

    std::tie(mesh, vpg) = gcs::makeHalfedgeAndGeometry(
        soup.polygons, soup.vertexCoordinates, true);
  }
  // ~MappingUtilityTest() {}
  // virtual void SetUp() {}
  // virtual void TearDown() {}

  std::unique_ptr<gcs::HalfedgeMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
};

TEST_F(MappingUtilityTest, MappingTest) {
  auto pos1 = EigenMap<double, 3>(vpg->inputVertexPositions);
  auto pos2 = EigenMap<double, 3>(vpg->inputVertexPositions.rawdata());

  ASSERT_EQ(vpg->inputVertexPositions.size(), pos1.rows());
  ASSERT_EQ(3, pos1.cols());
  ASSERT_EQ(vpg->inputVertexPositions.size(), pos2.rows());
  ASSERT_EQ(3, pos2.cols());

  // auto pos3 = EigenMap(vpg->inputVertexPositions);
  auto pos3 = vpg->inputVertexPositions.toMappedVector();
  auto pos4 = EigenMap(vpg->inputVertexPositions.rawdata());

  ASSERT_EQ(vpg->inputVertexPositions.size(), pos3.rows());
  ASSERT_EQ(1, pos3.cols());
  ASSERT_EQ(vpg->inputVertexPositions.size(), pos4.rows());
  ASSERT_EQ(1, pos4.cols());

  auto pos5 = FlattenedEigenMap<double, 3>(vpg->inputVertexPositions);
  auto pos6 = FlattenedEigenMap<double, 3>(vpg->inputVertexPositions.rawdata());

  ASSERT_EQ(3 * vpg->inputVertexPositions.size(), pos5.rows());
  ASSERT_EQ(1, pos5.cols());
  ASSERT_EQ(3 * vpg->inputVertexPositions.size(), pos6.rows());
  ASSERT_EQ(1, pos6.cols());

  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    std::size_t xidx = 3 * i;
    std::size_t yidx = 3 * i + 1;
    std::size_t zidx = 3 * i + 2;
    // check x values
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos1(i, 0));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos2(i, 0));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos3(i)[0]);
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos4(i)[0]);
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos5(xidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos6(xidx));

    // check y values
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos1(i, 1));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos2(i, 1));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos3(i)[1]);
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos4(i)[1]);
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos5(yidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos6(yidx));

    // check z values
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos1(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos2(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos3(i)[2]);
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos4(i)[2]);
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos5(zidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos6(zidx));

    vpg->inputVertexPositions[i][0] += 1;
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos1(i, 0));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos2(i, 0));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos3(i)[0]);
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos4(i)[0]);
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos5(xidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos6(xidx));

    pos1(i, 1) -= 3.14;
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos1(i, 1));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos2(i, 1));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos3(i)[1]);
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos4(i)[1]);
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos5(yidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos6(yidx));

    pos3(i)[2] *= 3.14;
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos1(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos2(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos3(i)[2]);
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos4(i)[2]);
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos5(zidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos6(zidx));

    pos5(zidx) /= 2.0;
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos1(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos2(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos3(i)[2]);
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos4(i)[2]);
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos5(zidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos6(zidx));
  }
}

TEST_F(MappingUtilityTest, MappingConstCorrectnessTest) {
  const gcs::VertexData<gc::Vector3> data(vpg->inputVertexPositions);
  const auto pos1 = EigenMap<double, 3>(data);
  const auto pos2 = EigenMap<double, 3>(data.rawdata());

  ASSERT_EQ(data.size(), pos1.rows());
  ASSERT_EQ(3, pos1.cols());
  ASSERT_EQ(data.size(), pos2.rows());
  ASSERT_EQ(3, pos2.cols());

  // auto pos3 = data.toMappedVector();
  const auto pos3 = EigenMap(data);
  const auto pos4 = EigenMap(data.rawdata());

  ASSERT_EQ(data.size(), pos3.rows());
  ASSERT_EQ(1, pos3.cols());
  ASSERT_EQ(data.size(), pos4.rows());
  ASSERT_EQ(1, pos4.cols());

  const auto pos5 = FlattenedEigenMap<double, 3>(data);
  const auto pos6 = FlattenedEigenMap<double, 3>(data.rawdata());

  ASSERT_EQ(3 * data.size(), pos5.rows());
  ASSERT_EQ(1, pos5.cols());
  ASSERT_EQ(3 * data.size(), pos6.rows());
  ASSERT_EQ(1, pos6.cols());
}
} // end namespace ddgsolver
