// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2020:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#include <iostream>

#include <gtest/gtest.h>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>

#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/util.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class MappingUtilityTest : public testing::Test {
protected:
  MappingUtilityTest() {
    std::vector<gc::Vector3> coords;
    std::vector<std::vector<std::size_t>> polygons;

    icosphere(coords, polygons, 2);

    gcs::PolygonSoupMesh soup(polygons, coords);
    soup.mergeIdenticalVertices();

    std::tie(mesh, vpg) =
        gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);
  }
  // ~MappingUtilityTest() {}
  // virtual void SetUp() {}
  // virtual void TearDown() {}

  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
};

TEST_F(MappingUtilityTest, MappingTest) {
  auto pos1 = EigenMap<double, 3>(vpg->inputVertexPositions);
  ASSERT_EQ(vpg->inputVertexPositions.size(), pos1.rows());
  ASSERT_EQ(3, pos1.cols());

  auto pos2 = FlattenedEigenMap<double, 3>(vpg->inputVertexPositions);
  ASSERT_EQ(3 * vpg->inputVertexPositions.size(), pos2.rows());
  ASSERT_EQ(1, pos2.cols());

  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    std::size_t xidx = 3 * i;
    std::size_t yidx = 3 * i + 1;
    std::size_t zidx = 3 * i + 2;
    // check x values
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos1(i, 0));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos2(xidx));

    // check y values
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos1(i, 1));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos2(yidx));

    // check z values
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos1(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos2(zidx));

    vpg->inputVertexPositions[i][0] += 1;
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos1(i, 0));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos2(xidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos1(i, 1));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos2(yidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos1(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos2(zidx));

    pos1(i, 1) -= 3.14;
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos1(i, 0));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos2(xidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos1(i, 1));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos2(yidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos1(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos2(zidx));

    pos2(zidx) /= 4;
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos1(i, 0));
    ASSERT_EQ(vpg->inputVertexPositions[i][0], pos2(xidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos1(i, 1));
    ASSERT_EQ(vpg->inputVertexPositions[i][1], pos2(yidx));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos1(i, 2));
    ASSERT_EQ(vpg->inputVertexPositions[i][2], pos2(zidx));
  }
}

TEST_F(MappingUtilityTest, MappingConstCorrectnessTest) {
  const gcs::VertexData<gc::Vector3> data(vpg->inputVertexPositions);
  const auto pos1 = EigenMap<double, 3>(data);
  ASSERT_EQ(data.size(), pos1.rows());
  ASSERT_EQ(3, pos1.cols());

  const auto pos2 = FlattenedEigenMap<double, 3>(data);
  ASSERT_EQ(3 * data.size(), pos2.rows());
  ASSERT_EQ(1, pos2.cols());
}
} // end namespace ddgsolver
