#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/util.h"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include <igl/cotmatrix.h>

int main() {
  Eigen::MatrixXd V(4, 2);
  V << 0, 0, 1, 0, 1, 1, 0, 1;
  Eigen::MatrixXi F(2, 3);
  F << 0, 1, 2, 0, 2, 3;
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  std::cout << "Hello, mesh: " << std::endl << L * V << std::endl;
  return 0;
}
