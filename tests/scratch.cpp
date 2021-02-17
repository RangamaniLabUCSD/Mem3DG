#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

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
  // Initialize unique ptr to mesh and geometry object
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrRefVpg;
  std::tie(ptrMesh, ptrVpg) = mem3dg::icosphere(1, 1);
  ptrRefVpg = ptrVpg->reinterpretTo(*ptrMesh);

  // List parameters
  double Kb = 8.22e-5, H0 = 0, sharpness = 10, Kst = 10, Ksl = 0, Kse = 0,
         epsilon = 15e-5, Bc = 40, gamma = 0, Vt = 0.7, Pam = 0, Kf = 0,
         conc = 25, height = 0, radius = 0.9, temp = 0, h = 5e-4, Kv = 5e-2,
         eta = 0, Ksg = 0.1;
  std::vector<double> pt = {1, 1, 1};
  std::vector<double> r_H0 = {100, 100};
  bool isProtein = false, isVertexShift = false, isReducedVolume = true,
       isLocalCurvature = false;
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);

  std::cout << "Initiating the system ...";
  mem3dg::Parameters p{Kb,    H0,  sharpness, r_H0, Ksg,    Kst,   Ksl, Kse,
                       Kv,    eta, epsilon,   Bc,   gamma,  Vt,    Pam, temp,
                       sigma, pt,  Kf,        conc, height, radius};
  mem3dg::System f(std::move(ptrMesh), std::move(ptrVpg), std::move(ptrRefVpg),
                   p, isReducedVolume, isProtein, isLocalCurvature,
                   isVertexShift);
  std::cout << "Finished!" << std::endl;

  std::cout << "Solving the system ..." << std::endl;
  double T = 3, eps = 0.002, closeZone = 1000, increment = 0, tSave = 1e-1,
         tMollify = 100;
  size_t verbosity = 0;

  mem3dg::ConjugateGradient integrator(f, h, true, T, tSave, eps, "./", "/traj.nc",
                           verbosity, true, 0.5, 1e-4, 0.01, false);
  integrator.integrate();

  return 0;
}