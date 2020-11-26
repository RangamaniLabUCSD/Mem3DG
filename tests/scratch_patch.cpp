#include <iostream>
#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

#include "mem3dg/solver/ddgsolver.h"
#include "mem3dg/solver/force.h"
#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/integrator.h"
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
  std::cout << "Loading input mesh " << inputMesh << " ...";
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;

  /*std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
  gcs::RichSurfaceMeshData::readMeshAndData(inputMesh); <- this returns no
  connectivity for UVsphere.ply ptrVpg = ptrRichData->getGeometry();*/
  std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(inputMesh);
  gcs::ManifoldSurfaceMesh &mesh = *ptrMesh;
  gcs::VertexPositionGeometry &vpg = *ptrVpg;
  std::cout << "Finished!" << std::endl;

  // To ensure that refVpg maps to mesh, take detour 
  // by first read and then construct
  std::string refMesh =
      "C://Users//Kieran//Dev//2020-Mem3DG-"
      "Applications//examples//patch_bud//input-file//patch.ply";
  std::cout << "Loading reference mesh " << refMesh << " ...";
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrRefMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrRefVpg_;
  std::tie(ptrRefMesh, ptrRefVpg_) = gcs::readManifoldSurfaceMesh(refMesh);
  gcs::VertexPositionGeometry *ptrRefVpg = new gcs::VertexPositionGeometry(
      mesh, gc::EigenMap<double, 3>(ptrRefVpg_->inputVertexPositions));

  std::cout << "Finished!" << std::endl;

  // size_t nSub = 0;
  // if (nSub > 0) {
  //   std::cout << "Subdivide input and reference mesh " << nSub
  //             << " time(s) ...";
  //   Eigen::Matrix<double, Eigen::Dynamic, 3> coords;
  //   Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> faces;

  //   Eigen::Matrix<double, Eigen::Dynamic, 3> refcoords;
  //   Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> reffaces;
  //   // ddgsolver::subdivide(ptrMesh, ptrVpg, nSub);
  //   // ddgsolver::subdivide(ptrRefMesh, ptrRefVpg, nSub);
  //   igl::loop(gc::EigenMap<double, 3>(ptrVpg->inputVertexPositions),
  //             ptrMesh->getFaceVertexMatrix<size_t>(), coords, faces, nSub);
  //   igl::loop(gc::EigenMap<double, 3>(ptrRefVpg->inputVertexPositions),
  //             ptrRefMesh->getFaceVertexMatrix<size_t>(), refcoords, reffaces,
  //             nSub);
  //   std::tie(ptrMesh, ptrVpg) =
  //       gcs::makeManifoldSurfaceMeshAndGeometry(coords, faces);
  //   std::tie(ptrRefMesh, ptrRefVpg) =
  //       gcs::makeManifoldSurfaceMeshAndGeometry(refcoords, reffaces);
  //   std::cout << "Finished!" << std::endl;
  // }

  gcs::RichSurfaceMeshData richData(*ptrMesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*ptrVpg);

  std::cout << "Initiating the system ...";
  /// physical parameters
  double Kb = 8.22e-5, H0 = 40, sharpness = 50, Kst = 0, Ksl = 0,
         Kse = 0, epsilon = 15e-5, Bc = 40, gamma = 3, Vt = 1, Kf = 0,
         conc = 25, height = 0, radius = 100, kt = 0, h = 1e-5, Kv = 0, eta = 0,
         Ksg = 0.05, mollifyFactor = 1e-3;

  std::vector<double> pt = {0, 0, 0};
  std::vector<double> r_H0 = {0.15, 0.15};
  bool isProtein = false, isTuftedLaplacian = false, isVertexShift = false;
  double sigma = sqrt(2 * gamma * kt / h);

  if (ptrMesh->hasBoundary() && (Vt != 1.0)) {
    Vt = 1.0;
    std::cout << "Geometry is a patch, so change Vt to 1.0!" << std::endl;
  }

  ddgsolver::Parameters p{Kb,  H0,    sharpness, r_H0,    Ksg,  Kst,    Ksl,
                          Kse, Kv,    eta,       epsilon, Bc,   gamma,  Vt,
                          kt,  sigma + 1e-18, pt,        Kf,      conc, height, radius};
  ddgsolver::Force f(*ptrMesh, *ptrVpg, *ptrRefVpg, richData, p, isProtein,
                     isTuftedLaplacian, mollifyFactor, isVertexShift);
  std::cout << "Finished!" << std::endl;

  std::cout << "Solving the system ..." << std::endl;
  double T = 3, eps = 0.002, closeZone = 1000, increment = 0, tSave = 1e-1,
         tMollify = 100, errorJumpLim = 600;
  std::string outputDir = "C://Users//Kieran//Desktop//";
  size_t verbosity = 2;
  ddgsolver::integration::velocityVerlet(f, h, T, eps, closeZone, increment, Kv,
                                         Ksg, tSave, tMollify, verbosity,
                                         inputMesh, outputDir, 0, errorJumpLim);

  delete ptrRefVpg;
  return 0;
}