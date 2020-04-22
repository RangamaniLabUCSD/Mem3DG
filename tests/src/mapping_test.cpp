

#include <iostream>

#include <gtest/gtest.h>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>

#include "ddgsolver/icosphere.h"
#include "ddgsolver/util.h"


namespace gc  = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class MappingUtilityTest : public testing::Test {};

TEST_F(MappingUtilityTest, ValidMappingTest) {
    std::vector<gc::Vector3> coords;
    std::vector<std::vector<std::size_t>> polygons;

    icosphere(coords, polygons, 2);

    gc::PolygonSoupMesh soup(polygons, coords);
    soup.mergeIdenticalVertices();

    std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
    std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
    std::tie(ptrmesh, ptrvpg) = gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

    auto& mesh = *ptrmesh;
    auto& vpg  = *ptrvpg;

    auto pos = mapVecToEigen<double, 3>(vpg.inputVertexPositions.rawdata());

    for (std::size_t i = 0; i < mesh.nVertices(); ++i){
        ASSERT_EQ(vpg.inputVertexPositions[i][0], pos(3*i));
        ASSERT_EQ(vpg.inputVertexPositions[i][1], pos(3*i+1));
        ASSERT_EQ(vpg.inputVertexPositions[i][2], pos(3*i+2));

        vpg.inputVertexPositions[i][0] += 1;
        ASSERT_EQ(vpg.inputVertexPositions[i][0], pos(3*i));

        pos(3*i) += 1;
        ASSERT_EQ(vpg.inputVertexPositions[i][0], pos(3*i));
    }
}