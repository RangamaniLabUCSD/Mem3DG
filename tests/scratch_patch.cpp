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
  // pybind11::scoped_interpreter guard{};
  std::string inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//"
                          "examples//patch_bud//input-"
                          "file//patch.ply";

  mem3dg::solver::Parameters p;
  p.bending.Kb = 8.22e-5;
  p.bending.Kbc = 8.22e-4;
  p.bending.H0c = 40;
  p.proteinDistribution.protein0.resize(2);
  p.proteinDistribution.protein0 << 0.15, 0.15;
  p.tension.Ksg = 0.05;
  p.tension.A_res = 0;
  p.osmotic.isPreferredVolume = true;
  p.osmotic.Kv = 0;
  p.adsorption.epsilon = 15e-5;
  p.proteinMobility = 40;
  p.point.pt.resize(3);
  p.point.pt << 0, 0, 0;
  p.variation.isProteinVariation = false;
  p.point.isFloatVertex = false;

  mem3dg::solver::System f(inputMesh, p, 0, false);

  double T = 3, h = 1e-5, eps = 0.002, closeZone = 1000, increment = 0,
         tSave = 1e-1, tMollify = 100, errorJumpLim = 600;
  std::string outputDir = "C://Users//Kieran//Desktop//";
  std::size_t verbosity = 2;

  mem3dg::solver::integrator::Euler integrator(f, h, T, tSave, eps, outputDir);
  integrator.isAdaptiveStep = true;
  integrator.trajFileName = "traj.nc";
  integrator.verbosity = verbosity;
  integrator.isBacktrack = true;
  integrator.rho = 0.5;
  integrator.c1 = 1e-4;
  integrator.integrate();

  return 0;
}
