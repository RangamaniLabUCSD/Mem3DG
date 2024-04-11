/*
 * Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
 *
 * Copyright 2024- The Mem3DG Authors
 * and the project initiators Cuncheng Zhu, Christopher T. Lee, and
 * Padmini Rangamani.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Please help us support Mem3DG development by citing the research
 * papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
 * for more information.
 */

#include <iostream>

#include <netcdf>

#include "mem3dg/mem3dg"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace nc = ::netCDF;

int main() {
  // std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  // std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  // std::tie(mesh, vpg) = mem3dg::icosphere(1, 0);

  // mem3dg::solver::MutableTrajFile f =
  // mem3dg::solver::MutableTrajFile::newFile(
  //     "/Users/ctlee/Mem3DG/build/bin/test.nc", true);

  // mem3dg::EigenVectorX3ur t1{mesh->getFaceVertexMatrix<std::uint32_t>()};
  // f.writeTopology(0, t1);
  // f.writeCoords(0, *vpg);
  // f.writeTime(0, 1.4);

  // std::tie(mesh, vpg) = mem3dg::icosphere(1, 1);
  // mem3dg::EigenVectorX3ur t2{mesh->getFaceVertexMatrix<std::uint32_t>()};
  // f.writeTopology(1, t2);
  // f.writeCoords(1, *vpg);
  // // f.writeVel(1, *vpg);

  // f.writeTime(1, 2.3);

  // f.close();

  // f.open("/Users/ctlee/Mem3DG/build/bin/test.nc", nc::NcFile::read);

  // for (std::size_t i = 0; i < 2; ++i) {
  //   std::cout << "Get Topology " << std::endl;
  //   std::cout << f.getTopology(i) << std::endl;
  //   std::cout << "Get coords " << std::endl;
  //   std::cout << f.getCoords(i) << std::endl;
  //   std::cout << f.getTime(i) << std::endl;
  // }

  // std::cout << "EOF" << std::endl;
  return 0;
}
