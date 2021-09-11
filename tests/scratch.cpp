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
  std::tie(mesh, vpg) = mem3dg::getCylinderMatrix(1, 16, 60, 7.5, 0);

  /// physical parameters
  p.proteinMobility = 0;
  p.temperature = 0;
  p.point.pt.resize(3);
  p.point.pt << 1, 0, 10;
  p.point.isFloatVertex = true;
  p.proteinDistribution.protein0.resize(1);
  p.proteinDistribution.protein0 << 1;
  p.boundary.shapeBoundaryCondition = "roller";
  p.boundary.proteinBoundaryCondition = "none";
  p.variation.isProteinVariation = false;
  p.variation.isShapeVariation = true;
  p.variation.radius = -1;
  p.bending.Kb = 0;
  p.bending.Kbc = 8.22e-5;
  p.bending.H0c = 0;
  p.tension.isConstantSurfaceTension = true;
  p.tension.Ksg = 1e-7;
  p.tension.A_res = 0;
  p.tension.At = -1;
  p.tension.lambdaSG = 0;
  p.adsorption.epsilon = 0;
  p.osmotic.isPreferredVolume = true;
  p.osmotic.isConstantOsmoticPressure = false;
  p.osmotic.Kv = 0.01;
  p.osmotic.V_res = 4 / 3 * 3.14;
  p.osmotic.n = 1;
  p.osmotic.Vt = 65.6201 * 0.5;
  p.osmotic.cam = -1;
  p.osmotic.lambdaV = 0;
  p.dirichlet.eta = 0;
  p.dpd.gamma = 0;
  p.external.Kf = 0.01;

  mem3dg::solver::MeshProcessor mP;
  mP.meshMutator.shiftVertex = false;
  mP.meshMutator.flipNonDelaunay = true;
  // mP.meshMutator.splitLarge = true;
  // mP.meshMutator.splitFat = true;
  mP.meshMutator.splitSkinnyDelaunay = true;
  mP.meshMutator.splitCurved = true;
  mP.meshMutator.curvTol = 0.01;
  mP.meshMutator.collapseSkinny = true;

  mem3dg::solver::System system(mesh, vpg, p, mP, 0);

  const double dt = 1, T = 1000000, eps = 1e-6, tSave = 5000, verbosity = 5;
  const std::string outputDir = "/tmp";

  mem3dg::solver::integrator::Euler integrator{system, dt,  T,
                                               tSave,  eps, outputDir};
  integrator.processMeshPeriod = 2000;
  integrator.isBacktrack = true;
  integrator.isAdaptiveStep = true;
  integrator.verbosity = verbosity;
  integrator.integrate();
  return 0;
}
