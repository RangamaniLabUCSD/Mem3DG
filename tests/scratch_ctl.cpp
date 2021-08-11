#include <iostream>

#include "mem3dg.h"
#include "mem3dg/mem3dg"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int main() {
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::tie(mesh, vpg) = mem3dg::icosphere(1, 0);

  mem3dg::solver::MutableTrajFile f = mem3dg::solver::MutableTrajFile::newFile(
      "/Users/ctlee/Mem3DG/build/bin/test.nc", true);

  mem3dg::EigenVectorX3ur t1 = mem3dg::getFaceVertexMatrix(*mesh);
  f.writeTopoFrame(0, t1);

  f.writeTime(0, 1.4);

  std::tie(mesh, vpg) = mem3dg::icosphere(1, 1);
  mem3dg::EigenVectorX3ur t2 = mem3dg::getFaceVertexMatrix(*mesh);
  f.writeTopoFrame(1, t2);

  f.writeTime(1, 2.3);

  auto top = f.readTopoFrame(0);
  std::cout << top << std::endl << std::endl;

  top = f.readTopoFrame(1);
  std::cout << top << std::endl;

  std::cout << "EOF" << std::endl;
  return 0;
}
