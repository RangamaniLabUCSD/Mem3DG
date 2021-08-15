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

#include "mem3dg/mem3dg"
#include "mem3dg/type_utilities.h"
#include <Eigen/Core>

#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace nc = ::netCDF;

class MutableTrajfileTest : public ::testing::Test {
public:
  MutableTrajfileTest() { f.createNewFile("test.nc", nc::NcFile::replace); }

  mem3dg::solver::MutableTrajFile f;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  mem3dg::EigenVectorX3ur t1, t2;
  mem3dg::EigenVectorX3dr g1, g2;
};

TEST_F(MutableTrajfileTest, ReadWriteTopologyCoordinates) {
  std::tie(mesh, vpg) = mem3dg::icosphere(1, 0);
  t1 = mesh->getFaceVertexMatrix<std::uint32_t>();
  g1 = gc::EigenMap<double, 3>(vpg->inputVertexPositions);

  f.writeTopology(0, t1);
  f.writeCoords(0, g1);

  std::tie(mesh, vpg) = mem3dg::icosphere(1, 1);
  f.writeTopology(1, *mesh);
  f.writeCoords(1, *vpg);

  t2 = mesh->getFaceVertexMatrix<std::uint32_t>();
  g2 = gc::EigenMap<double, 3>(vpg->inputVertexPositions);

  mem3dg::EigenVectorX3ur top;
  top = f.getTopology(0);
  ASSERT_EQ(top, t1);

  mem3dg::EigenVectorX3dr coords;
  coords = f.getCoords(0);
  ASSERT_EQ(coords, g1);

  top = f.getTopology(1);
  ASSERT_EQ(top, t2);

  coords = f.getCoords(1);
  ASSERT_EQ(coords, g2);
}

#endif
