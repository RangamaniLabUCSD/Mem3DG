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

  /// physical parameters
  double Kb = 8.22e-5, Kbc = 10 * 8.22e-5, H0 = 40, Kst = 0, Ksl = 0, Kse = 0,
         epsilon = 15e-5, Bc = 40, gamma = 3, Vt = 1, cam = 0, Kf = 0,
         conc = 25, height = 0, radius = 100, temp = 0, h = 1e-5, Kv = 0,
         eta = 0, Ksg = 0.05, A_res = 0, V_res = 0, n = 1;

  EigenVectorX1d pt(3), r_H0(2);
  pt << 0, 0, 0;
  r_H0 << 0.15, 0.15;

  mem3dg::solver::Parameters p;
  p.bending.Kb = Kb;
  p.bending.Kbc = Kbc;
  p.bending.H0c = H0;
  p.protein0 = r_H0;
  p.tension.Ksg = Ksg;
  p.tension.A_res = A_res;
  p.osmotic.Kv = Kv;
  p.osmotic.V_res = V_res;
  p.osmotic.Vt = Vt;
  p.osmotic.cam = cam;
  p.osmotic.n = n;
  p.osmotic.isPreferredVolume = true;
  p.dirichlet.eta = eta;
  p.adsorption.epsilon = epsilon;
  p.Bc = Bc;
  p.dpd.gamma = gamma;
  p.temp = temp;
  p.pt = pt;
  p.external.Kf = Kf;
  p.external.conc = conc;
  p.external.height = height;
  p.radius = radius;
  p.variation.isProteinVariation = false;

  mem3dg::solver::Options o;
  o.isFloatVertex = false;

  mem3dg::solver::System f(inputMesh, p, o, 0, false);

  double T = 3, eps = 0.002, closeZone = 1000, increment = 0, tSave = 1e-1,
         tMollify = 100, errorJumpLim = 600;
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
