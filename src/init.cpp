// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2020:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/utilities/vector3.h"
#include <stdexcept>
#ifdef MEM3DG_WITH_NETCDF
#include "mem3dg/solver/trajfile.h"
#endif
#include "mem3dg/solver/system.h"

#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include <polyscope/polyscope.h>

#include <csignal>
#include <iomanip>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {

#ifdef MEM3DG_WITH_NETCDF

void System::mapContinuationVariables(std::string trajFile, int startingFrame) {

  // Open netcdf file
  mem3dg::TrajFile fd = mem3dg::TrajFile::openReadOnly(trajFile);
  fd.getNcFrame(startingFrame);

  // Check consistent topology for continuation
  if (!(mesh->nFaces() == fd.getTopology().rows() &&
        mesh->nVertices() == fd.getCoords(startingFrame).rows())) {
    throw std::logic_error(
        "Topology for continuation parameters mapping is not consistent!");
  }

  // Map continuation variables
  time = fd.getTime(startingFrame);
  proteinDensity.raw() = fd.getProteinDensity(startingFrame);
  gc::EigenMap<double, 3>(vel) = fd.getVelocity(startingFrame);

  // Recompute cached variables
  updateVertexPositions();
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readTrajFile(std::string trajFile, int startingFrame, size_t nSub) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> referenceMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> referenceVpg;
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;

  mem3dg::TrajFile fd = mem3dg::TrajFile::openReadOnly(trajFile);
  fd.getNcFrame(startingFrame);
  std::tie(mesh, vpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      fd.getCoords(startingFrame), fd.getTopology());
  std::cout << "Loaded input mesh from " << trajFile << " of frame "
            << startingFrame << std::endl;

  /// Load reference geometry ptrRefVpg onto ptrMesh object
  std::tie(referenceMesh, referenceVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(fd.getRefcoordinate(),
                                              fd.getTopology());
  std::cout << "Loaded reference mesh" << std::endl;

  /// Subdivide the mesh and geometry objects
  if (nSub > 0) {
    mem3dg::loopSubdivide(mesh, vpg, nSub);
    mem3dg::loopSubdivide(referenceMesh, referenceVpg, nSub);
    std::cout << "Subdivided input and reference mesh " << nSub << " time(s)"
              << std::endl;
  }

  // reinterpret referenceVpg to mesh instead of referenceMesh
  refVpg = referenceVpg->reinterpretTo(*mesh);

  return std::make_tuple(std::move(mesh), std::move(vpg), std::move(refVpg));
}
#endif

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMeshes(std::string inputMesh, std::string refMesh, size_t nSub) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> referenceMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> referenceVpg;
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) = gcs::readManifoldSurfaceMesh(inputMesh);
  std::cout << "Loaded input mesh " << inputMesh << std::endl;

  // Load input reference mesh and geometry
  std::tie(referenceMesh, referenceVpg) = gcs::readManifoldSurfaceMesh(refMesh);
  std::cout << "Loaded reference mesh " << refMesh << std::endl;

  // Check consistent topology
  if (!(mesh->nVertices() == referenceMesh->nVertices() &&
        mesh->nEdges() == referenceMesh->nEdges() &&
        mesh->nFaces() == referenceMesh->nFaces())) {
    throw std::logic_error(
        "Topology of input mesh and reference mesh is not consistent! If not "
        "referencing a mesh, please have the option isRefMesh on and have "
        "input mesh as the duplicated argument!");
  }

  // Subdivide the mesh and geometry objects
  if (nSub > 0) {
    // mem3dg::subdivide(mesh, vpg, nSub);
    // mem3dg::subdivide(ptrRefMesh, ptrRefVpg, nSub);
    mem3dg::loopSubdivide(mesh, vpg, nSub);
    mem3dg::loopSubdivide(referenceMesh, referenceVpg, nSub);
    std::cout << "Subdivided input and reference mesh " << nSub << " time(s)"
              << std::endl;
  }

  // reinterpret referenceVpg to mesh instead of referenceMesh.
  refVpg = referenceVpg->reinterpretTo(*mesh);

  return std::make_tuple(std::move(mesh), std::move(vpg), std::move(refVpg));
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMeshes(Eigen::Matrix<double, Eigen::Dynamic, 3> topologyMatrix,
                   Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix,
                   Eigen::Matrix<double, Eigen::Dynamic, 3> refVertexMatrix,
                   size_t nSub) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> referenceMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> referenceVpg;
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(vertexMatrix, topologyMatrix);
  std::cout << "Loaded input mesh " << std::endl;

  // Load input reference mesh and geometry
  std::tie(referenceMesh, referenceVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(refVertexMatrix, topologyMatrix);
  std::cout << "Loaded reference mesh " << std::endl;

  // Subdivide the mesh and geometry objects
  if (nSub > 0) {
    // mem3dg::subdivide(mesh, vpg, nSub);
    // mem3dg::subdivide(ptrRefMesh, ptrRefVpg, nSub);
    mem3dg::loopSubdivide(mesh, vpg, nSub);
    mem3dg::loopSubdivide(referenceMesh, referenceVpg, nSub);
    std::cout << "Subdivided input and reference mesh " << nSub << " time(s)"
              << std::endl;
  }

  // reinterpret referenceVpg to mesh instead of referenceMesh.
  refVpg = referenceVpg->reinterpretTo(*mesh);

  return std::make_tuple(std::move(mesh), std::move(vpg), std::move(refVpg));
}

void System::checkParametersAndOptions() {

  if (mesh->hasBoundary()) {
    O.isOpenMesh = true;
  }

  // check validity of parameters / options
  if ((O.isEdgeFlip || O.isGrowMesh) && O.isRefMesh) {
    throw std::logic_error(
        "Topology changes are not compatible with reference mesh!");
  }

  if (P.Kst != 0 && !O.isRefMesh) {
    throw std::logic_error("For topology changing simulation, conformal mesh "
                           "regularization Kst cannot be applied!");
  }

  if (mesh->hasBoundary()) {
    if (O.isReducedVolume || P.Vt != -1.0 || P.cam != 0.0) {
      throw std::logic_error(
          "For open boundary simulation, isReducedVolume has to be false, Vt "
          "has to be -1, and cam has to be 0 ");
    }
  }

  if (!O.isLocalCurvature) {
    if (P.eta != 0) {
      throw std::logic_error(
          "line tension eta has to be 0 for nonlocal curvature!");
    }
    if (P.r_H0 != std::vector<double>({-1, -1})) {
      throw std::logic_error("r_H0 has to be {-1, -1} for nonlocal curvature!");
    }
    if (P.Kbc != -1) {
      throw std::logic_error("Kbc has to be set to -1 for nonlocal curvature!");
    }
  }

  if (O.isReducedVolume) {
    if (P.cam != -1) {
      throw std::logic_error("ambient concentration cam has to be -1 for "
                             "reduced volume parametrized simulation!");
    }
  } else {
    if (P.Vt != -1) {
      throw std::logic_error("reduced volume Vt has to be -1 for "
                             "ambient pressure parametrized simulation!");
    }
  }

  if (!O.isProtein) {
    if (P.epsilon != -1 || P.Bc != -1) {
      throw std::logic_error("Binding constant Bc and binding energy "
                             "epsilon has to be both -1 for "
                             "protein binding disabled simulation!");
    }
  } else {
    if (O.isLocalCurvature) {
      throw std::logic_error("Local curvature should be deactivated with "
                             "protein binding activated!");
    }
  }

  if (P.Kf == 0) {
    if (P.conc != -1 || P.height != 0) {
      throw std::logic_error("With no external force, its concentration "
                             "should be disabled (=-1) "
                             "and prescribed height should be set to 0!");
    }
  }

  if (P.pt.size() > 3) {
    throw std::logic_error(
        "Length of p.pt cannnot exceed 3! Instruction: (Length=1) => (vertex "
        "index); (Length=2) => ([x,y] coordinate); (Length=3) => ([x,y,z] "
        "coordinate)");
  }

  if (P.pt.size() == 1 && O.isFloatVertex) {
    throw std::logic_error(
        "To have Floating vertex, one must specify vertex by coordinate!");
  }
}

void System::pcg_test() {
  // Generate a normal distribution around that mean
  std::normal_distribution<> normal_dist(0, 2);

  // Make a copy of the RNG state to use later
  pcg32 rng_checkpoint = rng;

  // Produce histogram
  std::map<int, int> hist;
  for (int n = 0; n < 10000; ++n) {
    ++hist[std::round(normal_dist(rng))];
  }
  std::cout << "Normal distribution around " << 0 << ":\n";
  for (auto p : hist) {
    std::cout << std::fixed << std::setprecision(1) << std::setw(2) << p.first
              << ' ' << std::string(p.second / 30, '*') << '\n';
  }

  // Produce information about RNG usage
  std::cout << "Required " << (rng - rng_checkpoint) << " random numbers.\n";
}

void System::initConstants() {
  // Initialize random number generator
  pcg_extras::seed_seq_from<std::random_device> seed_source;
  rng = pcg32(seed_source);

  // define local vpg
  const auto &localVpg = O.isRefMesh ? refVpg : vpg;
  localVpg->requireEdgeLengths();
  localVpg->requireFaceAreas();

  // Initialize V-E distribution matrix for line tension calculation
  if (P.eta != 0) {
    D = localVpg->d0.transpose();
    for (int k = 0; k < D.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it) {
        it.valueRef() = 0.5;
      }
    }
  }

  // Find "the" vertex
  findTheVertex(*localVpg);

  // Initialize the constant mask based on distance from the point specified
  mask.raw() =
      (heatMethodDistance(*localVpg, thePoint.nearestVertex()).raw().array() <
       P.radius)
          .matrix();

  // Mask boundary element
  if (mesh->hasBoundary()) {
    boundaryMask(*mesh, mask.raw());
  }

  // Explicitly cached the reference face areas data
  refFaceAreas = localVpg->faceAreas;

  // Explicitly cached the reference edge length data
  refEdgeLengths = localVpg->edgeLengths;

  // Initialize the constant target surface (total mesh) area
  refSurfaceArea =
      O.isOpenMesh ? computeProjectedArea(*localVpg) : refFaceAreas.raw().sum();

  // Initialize the constant target mean face area
  if (!O.isRefMesh || O.isGrowMesh) {
    meanTargetFaceArea = refFaceAreas.raw().sum() / mesh->nFaces();
  }

  // Initialize the constant target mean edge length
  if (!O.isRefMesh) {
    meanTargetEdgeLength = refEdgeLengths.raw().sum() / mesh->nEdges();
  }

  // Initialize the target constant cross length ration
  if (O.isRefMesh) {
    targetLcrs = computeLengthCrossRatio(*localVpg);
  }

  // Initialize the constant reference volume
  refVolume = O.isOpenMesh ? 0.0
                           : std::pow(refSurfaceArea / constants::PI / 4, 1.5) *
                                 (4 * constants::PI / 3);

  // Initialize the constant spontaneous curvature
  H0.raw().setConstant(mesh->nVertices(), 1, P.H0);

  // Initialize the constant bending rigidity
  Kb.raw().setConstant(mesh->nVertices(), 1, P.Kb);
}

void System::updateVertexPositions() {

  // refresh cached quantities after regularization
  vpg->refreshQuantities();

  // EigenMap commonly used matrices
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg->vertexNormals);
  auto positions = gc::EigenMap<double, 3>(vpg->inputVertexPositions);

  // initialize/update Laplacian matrix
  M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();

  // initialize/update mean curvature 1. cotan laplacian definition, 2. dihedral
  // definition
  H.raw() =
      O.isLaplacianMeanCurvature
          ? rowwiseDotProduct(M_inv * L * positions / 2.0, vertexAngleNormal_e)
          : M_inv * vpg->vertexMeanCurvatures.raw();

  // initialize/update Gaussian curvature
  K.raw() = M_inv * vpg->vertexGaussianCurvatures.raw();

  /// initialize/update enclosed volume
  volume = getMeshVolume(*mesh, *vpg, gc::Vector3{0, 0, 0});

  // initialize/update total surface area
  surfaceArea = vpg->faceAreas.raw().sum();

  // update reference area by projecting to xy plane
  if (O.isOpenMesh) {
    refSurfaceArea = computeProjectedArea(*vpg);
  }

  // initialize/update spontaneous curvature (protein
  // binding)
  if (O.isProtein) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    H0.raw() =
        (P.H0 * proteinDensitySq.array() / (1 + proteinDensitySq.array()))
            .matrix();
  }

  // initialize/update line tension (on dual edge)
  if (P.eta != 0) {
    // scale the dH0 such that it is integrated over the edge
    // this is under the case where the resolution is low. This is where the
    // extra vpg->edgeLength comes from!!!
    // WIP The unit of line tension is in force*length (e.g. XXNewton)
    lineTension.raw() = P.eta * vpg->edgeLengths.raw().array() *
                        (vpg->d0 * H0.raw()).cwiseAbs().array();
    // lineTension.raw() = P.eta * (vpg->d0 * H0.raw()).cwiseAbs().array();
  }

  // initialize/update external force
  if (P.Kf != 0) {
    computeExternalPressure();
  }

  // initialize/update the vertex position of the last
  // iteration
  pastPositions = vpg->inputVertexPositions;
}

void System::findTheVertex(gcs::VertexPositionGeometry &vpg) {
  if (P.pt.size() == 1) {
    thePoint.vertex = mesh->vertex((std::size_t)P.pt[0]);
  } else if (P.pt.size() == 2) {
    auto closestVertex = closestVertexToPt(*mesh, vpg, P.pt);

    for (gcs::Halfedge he : closestVertex.outgoingHalfedges()) {
      gc::Vector2 v1{vpg.inputVertexPositions[he.vertex()].x,
                     vpg.inputVertexPositions[he.vertex()].y};
      gc::Vector2 v2{vpg.inputVertexPositions[he.next().vertex()].x,
                     vpg.inputVertexPositions[he.next().vertex()].y};
      gc::Vector2 v3{vpg.inputVertexPositions[he.next().next().vertex()].x,
                     vpg.inputVertexPositions[he.next().next().vertex()].y};
      gc::Vector2 v{P.pt[0], P.pt[1]};
      gc::Vector3 baryCoords_ = cartesianToBarycentric(v1, v2, v3, v);
      
      if (baryCoords_.x > 0 && baryCoords_.y > 0 && baryCoords_.z > 0) {
        thePoint = gcs::SurfacePoint(
            he.face(), correspondBarycentricCoordinates(baryCoords_, he));
        break;
      }

    }
  }
}

} // namespace mem3dg