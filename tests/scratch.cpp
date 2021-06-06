#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

#include "mem3dg/constants.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg.h"
#include "mem3dg/solver/mesh.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/type_utilities.h"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

using EigenVectorX1D = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EigenVectorX1D_i = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using EigenVectorX3D =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenTopVec =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

int main() {
  // pybind11::scoped_interpreter guard{};
  std::string inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//"
                          "results//bud//asymm//testTraj//frame360.ply";

  /// physical parameters
  mem3dg::Parameters p;
  p.Kb = 8.22e-5;
  p.Kbc = 8.22e-4;
  p.H0c = 6;
  p.protein0 = EigenVectorX1D(1);
  p.protein0 << 1;
  p.eta = 0;
  p.Ksg = 2e-2;
  p.A_res = 0;
  p.Kst = 0; // 2e-6;
  p.Ksl = 1e-7;
  p.Kse = 1e-7;
  p.epsilon = -1;
  p.Bc = -1;
  p.Kv = 1;
  p.V_res = 0;
  p.Vt = -1;
  p.cam = 0;
  p.Kf = 0;
  p.conc = -1;
  p.height = 0;
  p.radius = 100000;
  p.gamma = 0;
  p.temp = 0;
  p.pt = EigenVectorX1D(2);
  p.pt << 0, 0;

  mem3dg::Options o;
  o.isProteinVariation = false;
  o.isReducedVolume = false;
  o.isEdgeFlip = true;
  o.isSplitEdge = true;
  o.isCollapseEdge = true;
  o.isVertexShift = false;
  o.isFloatVertex = true;

  mem3dg::System f(inputMesh, p, o, 0, false);

  double h = 0.05, T = 4076, eps = 0, tSave = 10, rho = 0.99, c1 = 0.0001,
         verbosity = 3, restartNum = 5;
  bool isAdaptiveStep = true, isAugmentedLagrangian = false,
       isBacktrack = false;
  std::string outputDir = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//"
                          "results//bud//asymm//testTraj";

  mem3dg::Euler integrator(f, h, T, tSave, eps, outputDir, isAdaptiveStep,
                           "traj.nc", verbosity, isBacktrack, rho, c1);
  integrator.integrate();

  return 0;
}
