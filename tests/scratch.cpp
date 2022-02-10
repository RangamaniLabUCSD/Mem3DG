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

  // Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> mesh;
  // Eigen::Matrix<double, Eigen::Dynamic, 3> vpg;
  mem3dg::solver::Parameters p;
  // std::tie(mesh, vpg) = mem3dg::getCylinderMatrix(1, 16, 60, 7.5, 0);
  std::string inputMesh = "/home/cuzhu/Mem3DG/tests/frame9.ply";

  /// physical parameters
  p.proteinMobility = 10;
  p.temperature = 0;

  p.point.pt.resize(2);
  p.point.pt << 0, 0;
  p.point.isFloatVertex = false;

  p.proteinDistribution.profile = "none";
  p.proteinDistribution.protein0.resize(1);
  p.proteinDistribution.protein0 << -1;
  p.proteinDistribution.lambdaPhi = 1e-7;

  p.boundary.shapeBoundaryCondition = "fixed";
  p.boundary.proteinBoundaryCondition = "pin";

  p.variation.isProteinVariation = true;
  p.variation.isShapeVariation = true;
  p.variation.radius = -1;

  p.bending.Kb = 8.22e-5;
  p.bending.Kbc = 2 * 8.22e-5;
  p.bending.H0c = -60;

  p.tension.isConstantSurfaceTension = false;
  p.tension.Ksg = 1;
  p.tension.A_res = 0;
  p.tension.At = 3.40904;
  p.tension.lambdaSG = 0;

  p.adsorption.epsilon = -1e-4;

  p.aggregation.chi = -2e-1;

  p.osmotic.isPreferredVolume = false;
  p.osmotic.isConstantOsmoticPressure = true;
  p.osmotic.Kv = 0;
  p.osmotic.V_res = 0;
  p.osmotic.n = 1;
  p.osmotic.Vt = -1;
  p.osmotic.cam = -1;
  p.osmotic.lambdaV = 0;

  p.dirichlet.eta = 0.00005;

  p.selfAvoidance.d = 0.01;
  p.selfAvoidance.mu = 1e-6;
  p.selfAvoidance.p = 0.1;

  p.dpd.gamma = 0;
  p.external.Kf = 0;

  mem3dg::solver::MeshProcessor mP;
  mP.meshMutator.shiftVertex = true;
  mP.meshMutator.flipNonDelaunay = true;
  // mP.meshMutator.splitLarge = true;
  mP.meshMutator.splitFat = true;
  mP.meshMutator.splitSkinnyDelaunay = true;
  mP.meshMutator.splitCurved = true;
  mP.meshMutator.curvTol = 0.003;
  mP.meshMutator.collapseSkinny = true;

  // mem3dg::solver::System system(mesh, vpg, p, mP, 0);
  mem3dg::solver::System system(inputMesh, p, mP, 0, 0, true);

  const double dt = 0.01, T = 1000000, eps = 1e-4, tSave = 2, verbosity = 5;
  const std::string outputDir = "/tmp";

  mem3dg::solver::integrator::Euler integrator{system, dt,  T,
                                               tSave,  eps, outputDir};
  integrator.processMeshPeriod = 0.1;
  integrator.isBacktrack = true;
  integrator.isAdaptiveStep = true;
  integrator.verbosity = verbosity;
  integrator.integrate();
  return 0;
}
