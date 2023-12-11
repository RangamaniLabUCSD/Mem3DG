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
template <typename T>
using EigenVectorX1_T = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using EigenVectorX1i = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using EigenVectorX3dr =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3ur =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3u = Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3>;
using EigenVectorX3sr =
    Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3s = Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>;

int main() {
  // pybind11::scoped_interpreter guard{};
  std::string inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//"
                          "examples//patch_bud//input-"
                          "file//patch.ply";

  EigenVectorX3sr mesh;
  EigenVectorX3dr vpg;
  std::tie(mesh, vpg) = mem3dg::getFaceAndVertexMatrix(inputMesh);

  mem3dg::solver::Parameters p;
  p.bending.Kb = 8.22e-5;
  p.bending.Kbc = 8.22e-4;
  p.bending.H0c = 40;
  auto constantSurfaceTensionModel = [](double area) {
    double tension = 0.05;
    double energy = tension * area;
    return std::make_tuple(tension, energy);
  };
  p.tension.form = constantSurfaceTensionModel;
  p.adsorption.epsilon = 15e-5;
  p.proteinMobility = 40;

  std::size_t notableVertex_index =
      mem3dg::getVertexClosestToEmbeddedCoordinate(
          vpg, std::array<double, 3>{0, 0, 0},
          std::array<bool, 3>{true, true, false});
  Eigen::Matrix<bool, Eigen::Dynamic, 1> notableVertex =
      Eigen::Matrix<bool, Eigen::Dynamic, 1>::Constant(vpg.rows(), false);
  notableVertex[notableVertex_index] = true;
  p.variation.isProteinVariation = false;
  mem3dg::solver::Geometry geometry(inputMesh, notableVertex);
  mem3dg::solver::System system(geometry, p);
  system.initialize();

  double T = 3, h = 1e-5, eps = 0.002, closeZone = 1000, increment = 0,
         tSave = 1e-1, tMollify = 100, errorJumpLim = 600;
  std::string outputDir = "C://Users//Kieran//Desktop//";

  mem3dg::solver::integrator::Euler integrator(system, h, T, tSave, eps,
                                               outputDir);
  integrator.ifAdaptiveStep = true;
  integrator.trajFileName = "traj.nc";
  integrator.isBacktrack = true;
  integrator.rho = 0.5;
  integrator.c1 = 1e-4;
  integrator.integrate();

  return 0;
}
