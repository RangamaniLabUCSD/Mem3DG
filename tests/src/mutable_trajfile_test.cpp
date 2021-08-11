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
#include <Eigen/Core>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class MutableTrajfileTest : public ::testing::Test {
public:
  MutableTrajfileTest() {
    std::tie(mesh, vpg) = mem3dg::icosphere(1, 0);
    t1 = mem3dg::getFaceVertexMatrix(*mesh);

    std::tie(mesh, vpg) = mem3dg::icosphere(1, 1);
    t2 = mem3dg::getFaceVertexMatrix(*mesh);

    f.createNewFile("test.nc", NC_NETCDF4 | NC_CLOBBER | NC_DISKLESS);
    // mem3dg::solver::MutableTrajFile f =
    //     mem3dg::solver::MutableTrajFile::newDisklessFile("./test.nc");
  }

  mem3dg::solver::MutableTrajFile f;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  mem3dg::EigenVectorX3ur t1, t2;
};

TEST_F(MutableTrajfileTest, ReadWriteTopology) {
  f.writeTopoFrame(0, t1);
  f.writeTopoFrame(1, t2);

  mem3dg::EigenVectorX3ur top;
  top = f.readTopoFrame(0);
  ASSERT_EQ(top, t1);

  top = f.readTopoFrame(1);
  ASSERT_EQ(top, t2);
}
