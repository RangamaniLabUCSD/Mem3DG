#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

#include "mem3dg/solver/constants.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/mem3dg.h"
#include "mem3dg/solver/mesh.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/util.h"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int main() {
  // pybind11::scoped_interpreter guard{};
  std::string inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//"
                          "examples//patch_bud//input-"
                          "file//patch.ply";
  std::string refMesh =
      "C://Users//Kieran//Dev//2020-Mem3DG-"
      "Applications//examples//patch_bud//input-file//patch.ply";

  /// physical parameters
  double Kb = 8.22e-5, Kbc = 10 * 8.22e-5, H0 = 40, Kst = 0, Ksl = 0, Kse = 0,
         epsilon = 15e-5, Bc = 40, gamma = 3, Vt = 1, Pam = 0, Kf = 0,
         conc = 25, height = 0, radius = 100, temp = 0, h = 1e-5, Kv = 0,
         eta = 0, Ksg = 0.05;

  std::vector<double> pt = {0, 0, 0};
  std::vector<double> r_H0 = {0.15, 0.15};
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);

  mem3dg::Parameters p{Kb,    Kbc, H0,      r_H0, Ksg,    Kst,   Ksl, Kse,
                       Kv,    eta, epsilon, Bc,   gamma,  Vt,    Pam, temp,
                       sigma, pt,  Kf,      conc, height, radius};

   mem3dg::Options o;
  o.isProtein = false;
  o.isVertexShift = false;
  o.isReducedVolume = true;
  o.isLocalCurvature = true;
  o.isEdgeFlip = false;
  o.isGrowMesh = false;
  o.isRefMesh = true;
  o.isFloatVertex = false;
  o.isLaplacianMeanCurvature = false;

  mem3dg::System f(inputMesh, refMesh, 0, p, o);

  double T = 3, eps = 0.002, closeZone = 1000, increment = 0, tSave = 1e-1,
         tMollify = 100, errorJumpLim = 600;
  std::string outputDir = "C://Users//Kieran//Desktop//";
  size_t verbosity = 2;

  mem3dg::Euler integrator(f, h, true, T, tSave, eps, outputDir, "/traj.nc",
                           verbosity, true, 0.5, 1e-4);
  integrator.integrate();

  return 0;
}