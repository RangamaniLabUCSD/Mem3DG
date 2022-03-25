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

#include "mem3dg/solver/system.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {
namespace solver {
#ifdef MEM3DG_WITH_NETCDF
void System::mapContinuationVariables(std::string trajFile, int startingFrame) {

  // Open netcdf file
  // TrajFile fd = TrajFile::openReadOnly(trajFile);
  // fd.getNcFrame(startingFrame);
  MutableTrajFile fd = MutableTrajFile::openReadOnly(trajFile);
  fd.getNcFrame(startingFrame);

  // Check consistent topology for continuation
  if (!(mesh->nFaces() == fd.getTopology(startingFrame).rows() &&
        mesh->nVertices() == fd.getCoords(startingFrame).rows())) {
    throw std::logic_error(
        "Topology for continuation parameters mapping is not consistent!");
  } else {
    // Map continuation variables
    time = fd.getTime(startingFrame);
    energy.time = time;
    frame = startingFrame;
    toMatrix(velocity) = fd.getVelocity(startingFrame);
    // F.toMatrix(vel_protein) = fd.getProteinVelocity(startingFrame);
    if (parameters.proteinDistribution.protein0.rows() == 1 &&
        parameters.proteinDistribution.protein0[0] == -1) {
      proteinDensity.raw() = fd.getProteinDensity(startingFrame);
    } else {
      throw std::logic_error("proteinDensity.protein0 has to be disabled "
                             "(=[-1]) for continuing simulations!");
    }
  }
}
#endif

void System::mapContinuationVariables(std::string plyFile) {
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh_local;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData_local;

  // Open mesh file
  std::tie(ptrMesh_local, ptrRichData_local) =
      gcs::RichSurfaceMeshData::readMeshAndData(plyFile);

  // Check consistent topology for continuation
  if (mesh->nFaces() != ptrMesh_local->nFaces() ||
      mesh->nVertices() != ptrMesh_local->nVertices()) {
    throw std::logic_error(
        "Topology for continuation parameters mapping is not consistent!");
  } else {
    // Map continuation variables
    time = std::stod(sliceString(plyFile, "t", "_"));
    energy.time = time;
    frame = std::stod(sliceString(plyFile, "f", "_"));
    if (parameters.proteinDistribution.protein0.rows() == 1 &&
        parameters.proteinDistribution.protein0[0] == -1) {
      proteinDensity =
          ptrRichData_local->getVertexProperty<double>("protein_density")
              .reinterpretTo(*mesh);
      // vel_protein =
      //     ptrRichData_local->getVertexProperty<double>("protein_velocity")
      //         .reinterpretTo(*mesh);
      ;
    } else {
      throw std::logic_error("proteinDensity.protein0 has to be disabled "
                             "(=[-1]) for continuing simulations!");
    }
  }
}

#ifdef MEM3DG_WITH_NETCDF
std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readTrajFile(std::string trajFile, int startingFrame) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;

  // TrajFile fd = TrajFile::openReadOnly(trajFile);
  // fd.getNcFrame(startingFrame);
  MutableTrajFile fd = MutableTrajFile::openReadOnly(trajFile);
  fd.getNcFrame(startingFrame);
  std::tie(mesh, vpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      fd.getCoords(startingFrame), fd.getTopology(startingFrame));

  /// Subdivide the mesh and geometry objects
  // if (nSub > 0) {
  //   mem3dg::loopSubdivide(mesh, vpg, nSub);
  //   std::cout << "Subdivided input and reference mesh " << nSub << " time(s)"
  //             << std::endl;
  // }

  return std::make_tuple(std::move(mesh), std::move(vpg));
}
#endif

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMeshes(std::string inputMesh) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) = gcs::readManifoldSurfaceMesh(inputMesh);

  // // Subdivide the mesh and geometry objects
  // if (nSub > 0) {
  //   // mem3dg::subdivide(mesh, vpg, nSub);
  //   mem3dg::loopSubdivide(mesh, vpg, nSub);
  //   std::cout << "Subdivided input mesh " << nSub << " time(s)" << std::endl;
  // }

  return std::make_tuple(std::move(mesh), std::move(vpg));
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMeshes(
    Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &topologyMatrix,
    Eigen::Matrix<double, Eigen::Dynamic, 3> &vertexMatrix) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(vertexMatrix, topologyMatrix);

  // Subdivide the mesh and geometry objects
  // if (nSub > 0) {
  //   // mem3dg::subdivide(mesh, vpg, nSub);
  //   mem3dg::loopSubdivide(mesh, vpg, nSub);
  //   std::cout << "Subdivided input mesh " << nSub << " time(s)" << std::endl;
  // }

  return std::make_tuple(std::move(mesh), std::move(vpg));
}

void System::initialize(std::size_t nMutation) {
  initializeConstants();
  meshProcessor.summarizeStatus();
  if (!meshProcessor.isMeshMutate && nMutation != 0) {
    mem3dg_runtime_message("mesh mutator not activated!");
  } else {
    mutateMesh(nMutation);
  }
  updateConfigurations();
}

void System::checkConfiguration() {

  isOpenMesh = mesh->hasBoundary();
  parameters.checkParameters(isOpenMesh, mesh->nVertices());
  meshProcessor.summarizeStatus();
  if (meshProcessor.isMeshMutate && !parameters.variation.isShapeVariation) {
    mem3dg_runtime_error("Mesh mutation operation not allowed for non shape "
                         "variation simulation");
  }
  if (!isOpenMesh && mesh->genus() != 0) {
    mem3dg_runtime_error(
        "Do not support closed mesh with nonzero number of genus!")
  }
  if (meshProcessor.isMeshRegularize &&
      ((mesh->nVertices() != meshProcessor.meshRegularizer.nVertex ||
        mesh->nEdges() != meshProcessor.meshRegularizer.nEdge ||
        mesh->nFaces() != meshProcessor.meshRegularizer.nFace))) {
    mem3dg_runtime_error("For topologically different reference mesh, mesh "
                         "regularization cannot be applied!");
  }
  if (parameters.point.pt.rows() == 2 && !isOpenMesh) {
    mem3dg_runtime_message(
        "specifying x-y coordinate on closed surface may"
        "lead to ambiguity! Please check by visualizing it first!");
  }
  if (parameters.selfAvoidance.mu != 0) {
    for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
      gc::Vertex vi{mesh->vertex(i)};
      gc::VertexData<bool> neighborList(*mesh, false);
      meshProcessor.meshMutator.markVertices(neighborList, vi,
                                             parameters.selfAvoidance.n);
      for (std::size_t j = i + 1; j < mesh->nVertices(); ++j) {
        if (neighborList[j])
          continue;
        gc::Vertex vj{mesh->vertex(j)};
        gc::Vector3 r =
            vpg->inputVertexPositions[vj] - vpg->inputVertexPositions[vi];
        double distance = gc::norm(r);
        if (distance < parameters.selfAvoidance.d)
          mem3dg_runtime_error(
              "Input mesh violates the self avoidance constraint!");
      }
    }
  }
}

void System::initializeConstants() {
  // InitialiinitializeConstantsber generator
  pcg_extras::seed_seq_from<std::random_device> seed_source;
  rng = pcg32(seed_source);

  // // Initialize V-E distribution matrix for line tension calculation
  // if (P.dirichlet.eta != 0) {
  //   D = localVpg->d0.transpose().cwiseAbs() / 2;
  //   // for (int k = 0; k < D.outerSize(); ++k) {
  //   //   for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it;
  //   ++it)
  //   {
  //   //     it.valueRef() = 0.5;
  //   //   }
  //   // }
  // }

  // Find "the" vertex
  if (parameters.point.isFloatVertex) {
    findFloatCenter(*vpg, geodesicDistance, 1e18);
  } else {
    findVertexCenter(*vpg, geodesicDistance, 1e18);
  }

  // Initialize const geodesic distance
  updateGeodesicsDistance();

  // Initialize the constant mask based on distance from the point specified
  prescribeMasks();

  // Initialize protein density
  prescribeProteinDensity();

  // Mask boundary elementF
  if (mesh->hasBoundary()) {
    boundaryForceMask(*mesh, forces.forceMask,
                      parameters.boundary.shapeBoundaryCondition);
    boundaryProteinMask(*mesh, forces.proteinMask,
                        parameters.boundary.proteinBoundaryCondition);
  }

  // initialize/update total surface area
  surfaceArea = vpg->faceAreas.raw().sum() + parameters.tension.A_res;
  if (!ifMute)
    std::cout << "area_init = " << surfaceArea << std::endl;

  // Initialize the constant reference volume
  if (!ifMute)
    std::cout << "Characteristic volume wrt to At = "
              << (isOpenMesh
                      ? parameters.osmotic.V_res
                      : std::pow(parameters.tension.At / constants::PI / 4,
                                 1.5) *
                            (4 * constants::PI / 3))
              << std::endl;

  /// initialize/update enclosed volume
  volume = getMeshVolume(*mesh, *vpg, true) + parameters.osmotic.V_res;
  if (!ifMute)
    std::cout << "vol_init = " << volume << std::endl;
}

void System::updateConfigurations() {

  // refresh cached quantities after regularization
  vpg->refreshQuantities();

  // compute face gradient of protein density
  if (parameters.dirichlet.eta != 0) {
    computeGradient(proteinDensity, proteinDensityGradient);
  }

  // Update protein density dependent quantities
  if (parameters.bending.relation == "linear") {
    H0.raw() = proteinDensity.raw() * parameters.bending.H0c;
    Kb.raw() = parameters.bending.Kb +
               parameters.bending.Kbc * proteinDensity.raw().array();
    Kd.raw() = parameters.bending.Kd +
               parameters.bending.Kdc * proteinDensity.raw().array();
  } else if (parameters.bending.relation == "hill") {
    Eigen::Matrix<double, Eigen::Dynamic, 1> proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    H0.raw() = (parameters.bending.H0c * proteinDensitySq.array() /
                (1 + proteinDensitySq.array()))
                   .matrix();
    Kb.raw() = (parameters.bending.Kb + parameters.bending.Kbc *
                                            proteinDensitySq.array() /
                                            (1 + proteinDensitySq.array()))
                   .matrix();
    Kd.raw() = (parameters.bending.Kd + parameters.bending.Kdc *
                                            proteinDensitySq.array() /
                                            (1 + proteinDensitySq.array()))
                   .matrix();
  } else {
    mem3dg_runtime_error("updateVertexPosition: P.relation is invalid option!");
  }

  /// initialize/update enclosed volume
  volume = getMeshVolume(*mesh, *vpg, true) + parameters.osmotic.V_res;

  // update global osmotic pressure
  if (parameters.osmotic.isPreferredVolume) {
    forces.osmoticPressure =
        -(parameters.osmotic.Kv * (volume - parameters.osmotic.Vt) /
              parameters.osmotic.Vt / parameters.osmotic.Vt +
          parameters.osmotic.lambdaV);
  } else if (parameters.osmotic.isConstantOsmoticPressure) {
    forces.osmoticPressure = parameters.osmotic.Kv;
  } else {
    forces.osmoticPressure =
        mem3dg::constants::i * mem3dg::constants::R * parameters.temperature *
        (parameters.osmotic.n / volume - parameters.osmotic.cam);
  }

  // initialize/update total surface area
  surfaceArea = vpg->faceAreas.raw().sum() + parameters.tension.A_res;

  // update global surface tension
  forces.surfaceTension = parameters.tension.isConstantSurfaceTension
                              ? parameters.tension.Ksg
                              : parameters.tension.Ksg *
                                        (surfaceArea - parameters.tension.At) /
                                        parameters.tension.At +
                                    parameters.tension.lambdaSG;
}
} // namespace solver
} // namespace mem3dg
