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

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector3.h"
#include "mem3dg/constants.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/mutable_trajfile.h"
#include <cmath>
#include <stdexcept>
#include <vector>
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
  }

  // Map continuation variables
  time = fd.getTime(startingFrame);
  energy.time = time;
  proteinDensity.raw() = fd.getProteinDensity(startingFrame);
  toMatrix(velocity) = fd.getVelocity(startingFrame);
  // F.toMatrix(vel_protein) = fd.getProteinVelocity(startingFrame);
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readTrajFile(std::string trajFile, int startingFrame,
                     std::size_t nSub) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;

  // TrajFile fd = TrajFile::openReadOnly(trajFile);
  // fd.getNcFrame(startingFrame);
  MutableTrajFile fd = MutableTrajFile::openReadOnly(trajFile);
  fd.getNcFrame(startingFrame);
  std::tie(mesh, vpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      fd.getCoords(startingFrame), fd.getTopology(startingFrame));
  std::cout << "Loaded input mesh from " << trajFile << " of frame "
            << startingFrame << std::endl;

  /// Subdivide the mesh and geometry objects
  if (nSub > 0) {
    mem3dg::loopSubdivide(mesh, vpg, nSub);
    std::cout << "Subdivided input and reference mesh " << nSub << " time(s)"
              << std::endl;
  }

  return std::make_tuple(std::move(mesh), std::move(vpg));
}
#endif

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMeshes(std::string inputMesh, std::size_t nSub) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) = gcs::readManifoldSurfaceMesh(inputMesh);
  std::cout << "Loaded input mesh " << inputMesh << std::endl;

  // Subdivide the mesh and geometry objects
  if (nSub > 0) {
    // mem3dg::subdivide(mesh, vpg, nSub);
    mem3dg::loopSubdivide(mesh, vpg, nSub);
    std::cout << "Subdivided input mesh " << nSub << " time(s)" << std::endl;
  }

  return std::make_tuple(std::move(mesh), std::move(vpg));
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMeshes(
    Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &topologyMatrix,
    Eigen::Matrix<double, Eigen::Dynamic, 3> &vertexMatrix, std::size_t nSub) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(vertexMatrix, topologyMatrix);
  std::cout << "Loaded input mesh " << std::endl;

  // Subdivide the mesh and geometry objects
  if (nSub > 0) {
    // mem3dg::subdivide(mesh, vpg, nSub);
    mem3dg::loopSubdivide(mesh, vpg, nSub);
    std::cout << "Subdivided input mesh " << nSub << " time(s)" << std::endl;
  }

  return std::make_tuple(std::move(mesh), std::move(vpg));
}

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

void System::saveRichData(std::string PathToSave, bool isJustGeometry) {

  if (isJustGeometry) {
    gcs::writeSurfaceMesh(*mesh, *vpg, PathToSave);
  } else {
    gcs::RichSurfaceMeshData richData(*mesh);
    richData.addMeshConnectivity();
    richData.addGeometry(*vpg);

    // write protein distribution
    richData.addVertexProperty("protein_density", proteinDensity);
    richData.addVertexProperty("velocity", forces.ontoNormal(velocity));

    // write bool
    gcs::VertexData<double> msk(*mesh);
    msk.fromVector(toMatrix(forces.forceMask).rowwise().sum());
    richData.addVertexProperty("force_mask", msk);
    richData.addVertexProperty("protein_mask", forces.proteinMask);
    // gcs::VertexData<int> mutMkr(*mesh);
    // mutMkr.fromVector(mutationMarker.raw().cast<int>());
    // richData.addVertexProperty("smoothing_mask", mutMkr);
    gcs::VertexData<int> tkr(*mesh);
    tkr.fromVector(thePointTracker.raw().cast<int>());
    richData.addVertexProperty("the_point", tkr);

    // write geometry
    gcs::VertexData<double> meanCurv(*mesh);
    meanCurv.fromVector(vpg->vertexMeanCurvatures.raw().array() /
                        vpg->vertexDualAreas.raw().array());
    richData.addVertexProperty("mean_curvature", meanCurv);
    gcs::VertexData<double> gaussCurv(*mesh);
    gaussCurv.fromVector(vpg->vertexGaussianCurvatures.raw().array() /
                         vpg->vertexDualAreas.raw().array());
    richData.addVertexProperty("gauss_curvature", gaussCurv);
    richData.addVertexProperty("spon_curvature", H0);

    // write pressures
    richData.addVertexProperty("bending_force", forces.bendingForce);
    richData.addVertexProperty("deviatoric_force", forces.deviatoricForce);
    richData.addVertexProperty("capillary_force", forces.capillaryForce);
    richData.addVertexProperty("line_tension_force", forces.lineCapillaryForce);
    richData.addVertexProperty("osmotic_force", forces.osmoticForce);
    richData.addVertexProperty("adsorption_force", forces.adsorptionForce);
    richData.addVertexProperty("aggregation_force", forces.aggregationForce);
    richData.addVertexProperty("external_force", forces.externalForce);
    richData.addVertexProperty("avoidance_force", forces.selfAvoidanceForce);
    richData.addVertexProperty("physical_force", forces.mechanicalForce);

    // write chemical potential
    richData.addVertexProperty("diffusion_potential",
                               forces.diffusionPotential);
    richData.addVertexProperty("bending_potential", forces.bendingPotential);
    richData.addVertexProperty("deviatoric_potential",
                               forces.deviatoricPotential);
    richData.addVertexProperty("adsorption_potential",
                               forces.adsorptionPotential);
    richData.addVertexProperty("aggregation_potential",
                               forces.aggregationPotential);
    richData.addVertexProperty("chemical_potential", forces.chemicalPotential);

    richData.write(PathToSave);
  }
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

  // // Initialize V-E distribution matrix for line tension calculation
  // if (P.dirichlet.eta != 0) {
  //   D = localVpg->d0.transpose().cwiseAbs() / 2;
  //   // for (int k = 0; k < D.outerSize(); ++k) {
  //   //   for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it)
  //   {
  //   //     it.valueRef() = 0.5;
  //   //   }
  //   // }
  // }

  // Find "the" vertex
  findThePoint(*vpg, geodesicDistanceFromPtInd, 1e18);

  // Initialize const geodesic distance
  gcs::HeatMethodDistanceSolver heatSolver(*vpg);
  geodesicDistanceFromPtInd = heatSolver.computeDistance(thePoint);

  // Initialize the constant mask based on distance from the point specified
  if (parameters.variation.radius != -1) {
    if (parameters.variation.radius >
            geodesicDistanceFromPtInd.raw().maxCoeff() ||
        parameters.variation.radius <
            geodesicDistanceFromPtInd.raw().minCoeff()) {
      mem3dg_runtime_error("initConstants: either all vertices or none is "
                           "included within integration disk, "
                           "set radius = -1 to disable!");
    }
    for (gcs::Vertex v : mesh->vertices()) {
      forces.forceMask[v] =
          (geodesicDistanceFromPtInd[v] < parameters.variation.radius)
              ? gc::Vector3{1, 1, 1}
              : gc::Vector3{0, 0, 0};
      forces.proteinMask[v] =
          (geodesicDistanceFromPtInd[v] < parameters.variation.radius) ? 1 : 0;
    }
  }

  // Initialize protein density
  if (parameters.proteinDistribution.protein0.size() == 1) {
    proteinDensity.raw().setConstant(
        mesh->nVertices(), 1, parameters.proteinDistribution.protein0[0]);
  } else if (parameters.proteinDistribution.protein0.rows() ==
             proteinDensity.raw().rows()) {
    proteinDensity.raw() = parameters.proteinDistribution.protein0;
  } else if (parameters.proteinDistribution.protein0.size() == 4) {
    std::array<double, 2> r_heter{parameters.proteinDistribution.protein0[0],
                                  parameters.proteinDistribution.protein0[1]};
    vpg->requireVertexTangentBasis();
    if (parameters.proteinDistribution.profile == "gaussian") {
      gaussianDistribution(
          proteinDensity.raw(), geodesicDistanceFromPtInd.raw(),
          vpg->inputVertexPositions -
              vpg->inputVertexPositions[thePoint.nearestVertex()],
          vpg->vertexTangentBasis[thePoint.nearestVertex()], r_heter);
    } else if (parameters.proteinDistribution.profile == "tanh") {
      tanhDistribution(proteinDensity.raw(), geodesicDistanceFromPtInd.raw(),
                       vpg->inputVertexPositions -
                           vpg->inputVertexPositions[thePoint.nearestVertex()],
                       vpg->vertexTangentBasis[thePoint.nearestVertex()],
                       parameters.proteinDistribution.tanhSharpness, r_heter);
    }
    vpg->unrequireVertexTangentBasis();
    proteinDensity.raw() *= parameters.proteinDistribution.protein0[2] -
                            parameters.proteinDistribution.protein0[3];
    proteinDensity.raw().array() += parameters.proteinDistribution.protein0[3];
  }

  // Mask boundary elementF
  if (mesh->hasBoundary()) {
    boundaryForceMask(*mesh, forces.forceMask,
                      parameters.boundary.shapeBoundaryCondition);
    boundaryProteinMask(*mesh, forces.proteinMask,
                        parameters.boundary.proteinBoundaryCondition);
  }

  // initialize/update total surface area
  surfaceArea = vpg->faceAreas.raw().sum() + parameters.tension.A_res;
  std::cout << "area_init = " << surfaceArea << std::endl;

  // Initialize the constant reference volume
  std::cout << "Characteristic volume wrt to At = "
            << (isOpenMesh
                    ? parameters.osmotic.V_res
                    : std::pow(parameters.tension.At / constants::PI / 4, 1.5) *
                          (4 * constants::PI / 3))
            << std::endl;

  /// initialize/update enclosed volume
  volume = getMeshVolume(*mesh, *vpg, true) + parameters.osmotic.V_res;
  std::cout << "vol_init = " << volume << std::endl;
}

void System::updateConfigurations(bool isUpdateGeodesics) {

  // refresh cached quantities after regularization
  vpg->refreshQuantities();

  // recompute floating "the vertex"
  if (parameters.point.isFloatVertex && isUpdateGeodesics) {
    findThePoint(
        *vpg, geodesicDistanceFromPtInd,
        3 * vpg->edgeLength(thePoint.nearestVertex().halfedge().edge()));
  }

  // update geodesic distance
  if (isUpdateGeodesics) {
    gcs::HeatMethodDistanceSolver heatSolver(*vpg);
    geodesicDistanceFromPtInd = heatSolver.computeDistance(thePoint);
  }

  // initialize/update external force
  if (parameters.external.Kf != 0 && isUpdateGeodesics) {
    prescribeExternalForce();
  }

  // update protein density
  if (parameters.proteinDistribution.protein0.rows() == 4 &&
      !parameters.variation.isProteinVariation && isUpdateGeodesics) {
    std::array<double, 2> r_heter{parameters.proteinDistribution.protein0[0],
                                  parameters.proteinDistribution.protein0[1]};
    vpg->requireVertexTangentBasis();
    if (parameters.proteinDistribution.profile == "gaussian") {
      gaussianDistribution(
          proteinDensity.raw(), geodesicDistanceFromPtInd.raw(),
          vpg->inputVertexPositions -
              vpg->inputVertexPositions[thePoint.nearestVertex()],
          vpg->vertexTangentBasis[thePoint.nearestVertex()], r_heter);
    } else if (parameters.proteinDistribution.profile == "tanh") {
      tanhDistribution(proteinDensity.raw(), geodesicDistanceFromPtInd.raw(),
                       vpg->inputVertexPositions -
                           vpg->inputVertexPositions[thePoint.nearestVertex()],
                       vpg->vertexTangentBasis[thePoint.nearestVertex()],
                       parameters.proteinDistribution.tanhSharpness, r_heter);
    }
    vpg->unrequireVertexTangentBasis();
    proteinDensity.raw() *= parameters.proteinDistribution.protein0[2] -
                            parameters.proteinDistribution.protein0[3];
    proteinDensity.raw().array() += parameters.proteinDistribution.protein0[3];
  }

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

  // initialize/update line tension (on dual edge)
  if (parameters.dirichlet.eta != 0 && false) {
    mem3dg_runtime_error(
        "updateVertexPosition: out of data implementation on line tension, "
        "shouldn't be called!");
    // scale the dH0 such that it is integrated over the edge
    // this is under the case where the resolution is low. This is where the
    // extra vpg->edgeLength comes from!!!
    // WIP The unit of line tension is in force*length (e.g. XXNewton)
    // F.lineTension.raw() = P.dirichlet.eta * vpg->edgeLengths.raw().array() *
    //                       (vpg->d0 * H0.raw()).cwiseAbs().array();
    // lineTension.raw() = P.dirichlet.eta * (vpg->d0 *
    // H0.raw()).cwiseAbs().array();
  }
}

double System::inferTargetSurfaceArea() {
  double targetArea;
  if (isOpenMesh) {
    targetArea = parameters.tension.A_res;
    for (gcs::BoundaryLoop bl : mesh->boundaryLoops()) {
      targetArea += computePolygonArea(bl, vpg->inputVertexPositions);
    }
  } else {
    targetArea = vpg->faceAreas.raw().sum();
  }
  return targetArea;
}

void System::findThePoint(gcs::VertexPositionGeometry &vpg,
                          gcs::VertexData<double> &geodesicDistance,
                          double range) {
  bool isUpdated = false;
  if (parameters.point.isFloatVertex) {
    switch (parameters.point.pt.rows()) {
    case 1: {
      throw std::logic_error(
          "To have Floating vertex, one must specify vertex by coordinate!");
      break;
    }
    case 2: {
      // Find the cloest vertex to the point in the x-y plane
      gcs::Vertex closestVertex = closestVertexToPt(
          *mesh, vpg, parameters.point.pt, geodesicDistance, range);
      double shortestDistance = 1e18;
      // loop over every faces around the vertex
      for (gcs::Halfedge he : closestVertex.outgoingHalfedges()) {
        if (he.isInterior()) {
          // specify vertex coordinates and the target coordinate on the face
          gc::Vector2 v1{vpg.inputVertexPositions[he.vertex()].x,
                         vpg.inputVertexPositions[he.vertex()].y};
          gc::Vector2 v2{vpg.inputVertexPositions[he.next().vertex()].x,
                         vpg.inputVertexPositions[he.next().vertex()].y};
          gc::Vector2 v3{vpg.inputVertexPositions[he.next().next().vertex()].x,
                         vpg.inputVertexPositions[he.next().next().vertex()].y};
          gc::Vector2 v{parameters.point.pt[0], parameters.point.pt[1]};
          // find the inverse barycentric mapping based on the cartesian vertex
          // coordinates
          gc::Vector3 baryCoords_ = cartesianToBarycentric(v1, v2, v3, v);

          if (baryCoords_.x > 0 && baryCoords_.y > 0 &&
              baryCoords_.z > 0) { // A. set the surface point when the point
                                   // lays within the triangle
            thePoint = gcs::SurfacePoint(
                he.face(), correspondBarycentricCoordinates(baryCoords_, he));
            isUpdated = true;
            break;
          } else { // B. avoid the floating point comparision, find the best by
                   // looping over the whole fan
            baryCoords_ =
                gc::componentwiseMax(baryCoords_, gc::Vector3{0, 0, 0});
            baryCoords_ /= gc::sum(baryCoords_);
            gcs::SurfacePoint someSurfacePoint(
                he.face(), correspondBarycentricCoordinates(baryCoords_, he));
            double distance =
                (gc::Vector2{parameters.point.pt[0], parameters.point.pt[1]} -
                 gc::Vector2{
                     someSurfacePoint.interpolate(vpg.inputVertexPositions).x,
                     someSurfacePoint.interpolate(vpg.inputVertexPositions).y})
                    .norm();
            if (distance < shortestDistance) {
              thePoint = someSurfacePoint;
              shortestDistance = distance;
              isUpdated = true;
            }
          }
        }
      }
      break;
    }
    case 3: {
      // initialize embedded point and the closest vertex
      gc::Vector3 embeddedPoint{parameters.point.pt[0], parameters.point.pt[1],
                                parameters.point.pt[2]};
      gcs::Vertex closestVertex = closestVertexToPt(
          *mesh, vpg, parameters.point.pt, geodesicDistance, range);
      gc::Vector3 vertexToPoint =
          embeddedPoint - vpg.inputVertexPositions[closestVertex];
      // initialize the surface point as the closest vertex
      // thePoint = gc::SurfacePoint(closestVertex);
      // double shortestDistance = vertexToPoint.norm();
      double shortestDistance = 1e18;
      // loop over every faces around the vertex
      for (gcs::Halfedge he : closestVertex.outgoingHalfedges()) {
        if (he.isInterior()) {
          // project the embedded point onto the face
          auto faceNormal = vpg.faceNormal(he.face());
          gc::Vector3 projectedEmbeddedPoint =
              embeddedPoint - gc::dot(vertexToPoint, faceNormal) * faceNormal;
          // determine the choice of coordinates used for inverse
          // barycentric mapping based on orientation of the face
          gc::Vector2 v1, v2, v3, v;
          if (abs(faceNormal.z) > std::sqrt(3) / 3) {
            v1 = gc::Vector2{vpg.inputVertexPositions[he.vertex()].x,
                             vpg.inputVertexPositions[he.vertex()].y};
            v2 = gc::Vector2{vpg.inputVertexPositions[he.next().vertex()].x,
                             vpg.inputVertexPositions[he.next().vertex()].y};
            v3 = gc::Vector2{
                vpg.inputVertexPositions[he.next().next().vertex()].x,
                vpg.inputVertexPositions[he.next().next().vertex()].y};
            v = gc::Vector2{projectedEmbeddedPoint.x, projectedEmbeddedPoint.y};
          } else if (abs(faceNormal.x) > std::sqrt(3) / 3) {
            v1 = gc::Vector2{vpg.inputVertexPositions[he.vertex()].y,
                             vpg.inputVertexPositions[he.vertex()].z};
            v2 = gc::Vector2{vpg.inputVertexPositions[he.next().vertex()].y,
                             vpg.inputVertexPositions[he.next().vertex()].z};
            v3 = gc::Vector2{
                vpg.inputVertexPositions[he.next().next().vertex()].y,
                vpg.inputVertexPositions[he.next().next().vertex()].z};
            v = gc::Vector2{projectedEmbeddedPoint.y, projectedEmbeddedPoint.z};
          } else {
            v1 = gc::Vector2{vpg.inputVertexPositions[he.vertex()].z,
                             vpg.inputVertexPositions[he.vertex()].x};
            v2 = gc::Vector2{vpg.inputVertexPositions[he.next().vertex()].z,
                             vpg.inputVertexPositions[he.next().vertex()].x};
            v3 = gc::Vector2{
                vpg.inputVertexPositions[he.next().next().vertex()].z,
                vpg.inputVertexPositions[he.next().next().vertex()].x};
            v = gc::Vector2{projectedEmbeddedPoint.z, projectedEmbeddedPoint.x};
          }
          // find the inverse barycentric mapping based on the cartesian
          // vertex coordinates
          gc::Vector3 baryCoords_ = cartesianToBarycentric(v1, v2, v3, v);
          // since might not find the perfect reflecting face, best we could
          // do within each triangle
          baryCoords_ = gc::componentwiseMax(baryCoords_, gc::Vector3{0, 0, 0});
          baryCoords_ /= gc::sum(baryCoords_);
          gcs::SurfacePoint someSurfacePoint(
              he.face(), correspondBarycentricCoordinates(baryCoords_, he));
          // compute optimum distance and set surface point
          double distance = (embeddedPoint - someSurfacePoint.interpolate(
                                                 vpg.inputVertexPositions))
                                .norm();
          if (distance < shortestDistance) {
            thePoint = someSurfacePoint;
            shortestDistance = distance;
            isUpdated = true;
          }
        }
      }
      break;
    }
    }
    // mark three vertex on the face
    thePointTracker.fill(false);
    thePointTracker[thePoint.face.halfedge().vertex()] = true;
    thePointTracker[thePoint.face.halfedge().next().vertex()] = true;
    thePointTracker[thePoint.face.halfedge().next().next().vertex()] = true;
  } else {
    switch (parameters.point.pt.rows()) {
    case 1: {
      // Assign surface point as the indexed vertex
      thePoint =
          gc::SurfacePoint(mesh->vertex((std::size_t)parameters.point.pt[0]));
      isUpdated = true;
      break;
    }
    case 2: {
      // Find the cloest vertex to the point in the x-y plane
      thePoint = gc::SurfacePoint(closestVertexToPt(
          *mesh, vpg, parameters.point.pt, geodesicDistance, range));
      isUpdated = true;
      break;
    }
    case 3: {
      thePoint = gc::SurfacePoint(closestVertexToPt(
          *mesh, vpg, parameters.point.pt, geodesicDistance, range));
      isUpdated = true;
      break;
    }
    }
    thePointTracker[thePoint.vertex] = true;
  }

  if (!isUpdated) {
    mem3dg_runtime_error("Surface point is not updated!");
  }
}
} // namespace solver
} // namespace mem3dg
