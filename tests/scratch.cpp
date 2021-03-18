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
#include "mem3dg/solver/visualization.h"

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

  // std::string inputMesh =
  // "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/bud.ply"; std::string
  // refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply";
  // std::string outputDir =
  // "/home/cuzhu/2020-Mem3DG-Applications/results/bud/coarse/Ksg5e-4_H5";
  std::string inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//"
                          "run//input-file//bud.ply";
  std::string refMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run/"
                        "/input-file//patch.ply";
  std::string outputDir = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//"
                          "results//bud//asymm//testTraj";
  // std::string inputMesh =
  // "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//newBUD.ply";
  // std::string refMesh =
  // "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//newBUD.ply";

  // Initialize unique ptr to mesh and geometry object
  // std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  // std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  // std::unique_ptr<gcs::VertexPositionGeometry> ptrRefVpg;
  // std::tie(ptrMesh, ptrVpg) = mem3dg::icosphere(1, 1);
  // ptrRefVpg = ptrVpg->reinterpretTo(*ptrMesh);

  // List parameters
  double Kb = 8.22e-5, Kbc = 10 * 8.22e-5, H0 = 4, Kst = 0, Ksl = 0, Kse = 0,
         epsilon = -1, Bc = -1, gamma = 0, Vt = 0.7, Pam = 0, Kf = 0, conc = 25,
         height = 0, radius = 1000, temp = 0, h = 5e-4, Kv = 5e-2, eta = 0,
         Ksg = 0.1;
  std::vector<double> pt = {1, 1, 1};
  std::vector<double> r_H0 = {0.5, 0.5};
  double sigma = sqrt(2 * gamma * mem3dg::constants::kBoltzmann * temp / h);

  std::cout << "Initiating the system ...";
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

  gcs::RichSurfaceMeshData richData(*f.mesh);
  // richData.addMeshConnectivity();
  // richData.addGeometry(*f.vpg);
  // std::string name{"//newBUD.ply"};
  // richData.write(outputDir + name);
  visualize(f);

  std::cout << "Finished!" << std::endl;

  //   std::cout << "Solving the system ..." << std::endl;
  //   double T = 3, eps = 0.002, closeZone = 1000, increment = 0, tSave = 1e-1,
  //          tMollify = 100;
  //   size_t verbosity = 0;

  //   mem3dg::ConjugateGradient integrator(f, h, true, T, tSave, eps, "./",
  //                                        "/traj.nc", verbosity, true, 0.5,
  //                                        1e-4, 0.01, false);
  //   integrator.integrate();

  return 0;
}