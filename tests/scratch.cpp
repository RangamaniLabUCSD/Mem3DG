#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

#include "mem3dg/mem3dg"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

using EigenVectorX1d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EigenVectorX1i = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using EigenVectorX3dr =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3ur =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

int main() {

  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> mesh;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vpg;
  mem3dg::solver::Parameters p;
  std::tie(mesh, vpg) = mem3dg::getIcosphereMatrix(1, 3);

  /// physical parameters
  p.bending.Kbc = 8.22e-5;
  p.tension.Ksg = 0.1;
  p.osmotic.isPreferredVolume = true;
  p.osmotic.Kv = 0.01;
  p.osmotic.Vt = 4.0 / 3.0 * mem3dg::constants::PI * 0.7;

  const double dt = 0.5, T = 50, eps = 0, tSave = 10, verbosity = 1;
  const std::string outputDir = "/tmp";

  mem3dg::solver::System f(mesh, vpg, p, 0);
  mem3dg::solver::integrator::Euler integrator{f, dt, T, tSave, eps, outputDir};
  integrator.trajFileName = "traj.nc";
  integrator.verbosity = verbosity;
  integrator.integrate();
  std::cout << "in main()" << std::endl;
  return 0;
}
