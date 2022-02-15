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

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"
#include <Eigen/Core>
#include <math.h>
#include <pcg_random.hpp>
#include <stdexcept>

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

gc::Vector3 System::cornerAngleGradient(gcs::Corner c, gcs::Vertex v) {
  gcs::Halfedge he = c.halfedge();
  gc::Vector3 n = vpg->faceNormals[c.face()];
  gc::Vector3 ej = vecFromHalfedge(he, *vpg);
  gc::Vector3 ei = vecFromHalfedge(he.next(), *vpg);
  gc::Vector3 ek = vecFromHalfedge(he.next().next(), *vpg);
  if (c.vertex() == v) { // vi
    gc::Vector3 grad_anglek = -gc::cross(n, ej).normalize() / gc::norm(ej);
    gc::Vector3 grad_anglej = -gc::cross(n, ek).normalize() / gc::norm(ek);
    return -(grad_anglek + grad_anglej);
  } else if (he.next().vertex() == v) { // vk
    return -gc::cross(n, ej).normalize() / gc::norm(ej);
  } else if (he.next().next().vertex() == v) { // vj
    return -gc::cross(n, ek).normalize() / gc::norm(ek);
  } else {
    mem3dg_runtime_error("Unexpected combination of corner and vertex!");
    return gc::Vector3{0, 0, 0};
  }
}

gc::Vector3 System::dihedralAngleGradient(gcs::Halfedge he, gcs::Vertex v) {
  double l = vpg->edgeLengths[he.edge()];
  if (he.edge().isBoundary()) {
    return gc::Vector3{0, 0, 0};
  } else if (he.vertex() == v) {
    return (vpg->halfedgeCotanWeights[he.next().next()] *
                vpg->faceNormals[he.face()] +
            vpg->halfedgeCotanWeights[he.twin().next()] *
                vpg->faceNormals[he.twin().face()]) /
           l;
  } else if (he.next().vertex() == v) {
    return (vpg->halfedgeCotanWeights[he.twin().next().next()] *
                vpg->faceNormals[he.twin().face()] +
            vpg->halfedgeCotanWeights[he.next()] *
                vpg->faceNormals[he.face()]) /
           l;
  } else if (he.next().next().vertex() == v) {
    return (-(vpg->halfedgeCotanWeights[he.next().next()] +
              vpg->halfedgeCotanWeights[he.next()]) *
            vpg->faceNormals[he.face()]) /
           l;
  } else {
    mem3dg_runtime_error("Unexpected combination of halfedge and vertex!");
    return gc::Vector3{0, 0, 0};
  }
}

std::tuple<gc::Vector3, gc::Vector3>
System::computeHalfedgeSchlafliVector(gcs::VertexPositionGeometry &vpg,
                                      gc::Halfedge &he) {
  std::size_t fID = he.face().getIndex();
  std::size_t heID_twin = he.twin().getIndex();
  std::size_t fID_he_twin = he.twin().face().getIndex();
  std::size_t heID_twin_next = he.twin().next().getIndex();
  std::size_t heID_he_next_next = he.next().next().getIndex();
  gc::Vertex vj = he.tipVertex();
  bool boundaryVertex = he.vertex().isBoundary();
  bool boundaryEdge = he.edge().isBoundary();
  bool interiorHalfedge = he.isInterior();
  bool interiorTwinHalfedge = he.twin().isInterior();
  gc::Vector3 schlafliVec1{0, 0, 0};
  gc::Vector3 schlafliVec2{0, 0, 0};
  if (!boundaryEdge) {
    schlafliVec1 =
        vpg.halfedgeCotanWeights[heID_he_next_next] * vpg.faceNormals[fID] +
        vpg.halfedgeCotanWeights[heID_twin_next] * vpg.faceNormals[fID_he_twin];
  }
  if (boundaryVertex && boundaryEdge) {
    schlafliVec2 = interiorHalfedge
                       ? (-(vpg.halfedgeCotanWeights[he] +
                            vpg.halfedgeCotanWeights[heID_he_next_next]) *
                          vpg.faceNormals[fID])
                       : (-(vpg.halfedgeCotanWeights[heID_twin] +
                            vpg.halfedgeCotanWeights[heID_twin_next]) *
                          vpg.faceNormals[fID_he_twin]);
  } else if (!boundaryVertex && vj.isBoundary()) {
    schlafliVec2 =
        vpg.halfedgeCotanWeights[heID_he_next_next] * vpg.faceNormals[fID] +
        vpg.halfedgeCotanWeights[heID_twin_next] * vpg.faceNormals[fID_he_twin];

    if (!he.next().edge().isBoundary())
      schlafliVec2 -= (vpg.halfedgeCotanWeights[he] +
                       vpg.halfedgeCotanWeights[heID_he_next_next]) *
                      vpg.faceNormals[fID];

    if (!he.twin().next().next().edge().isBoundary())
      schlafliVec2 -= (vpg.halfedgeCotanWeights[heID_twin] +
                       vpg.halfedgeCotanWeights[heID_twin_next]) *
                      vpg.faceNormals[fID_he_twin];
  } else {
    schlafliVec2 =
        -(vpg.halfedgeCotanWeights[he] * vpg.faceNormals[fID] +
          vpg.halfedgeCotanWeights[heID_twin] * vpg.faceNormals[fID_he_twin]);
  }
  return std::make_tuple(schlafliVec1, schlafliVec2);
}

gc::Vector3
System::computeHalfedgeGaussianCurvatureVector(gcs::VertexPositionGeometry &vpg,
                                               gc::Halfedge &he) {
  gc::Vector3 gaussVec{0, 0, 0};
  if (!he.edge().isBoundary()) {
    // gc::Vector3 eji{} = -vecFromHalfedge(he, *vpg);
    gaussVec = 0.5 * vpg.edgeDihedralAngles[he.edge()] *
               (-vecFromHalfedge(he, vpg)).unit();
  }
  return gaussVec;
}

gc::Vector3
System::computeHalfedgeMeanCurvatureVector(gcs::VertexPositionGeometry &vpg,
                                           gc::Halfedge &he) {
  std::size_t fID = he.face().getIndex();
  std::size_t fID_he_twin = he.twin().face().getIndex();
  bool interiorHalfedge = he.isInterior();
  bool interiorTwinHalfedge = he.twin().isInterior();
  gc::Vector3 areaGrad{0, 0, 0};
  if (interiorHalfedge) {
    areaGrad +=
        0.25 * gc::cross(vpg.faceNormals[fID], vecFromHalfedge(he.next(), vpg));
  }
  if (interiorTwinHalfedge) {
    areaGrad += 0.25 * gc::cross(vpg.faceNormals[fID_he_twin],
                                 vecFromHalfedge(he.twin().next().next(), vpg));
  }
  return areaGrad / 2;
}

gc::Vector3
System::computeHalfedgeVolumeVariationVector(gcs::VertexPositionGeometry &vpg,
                                             gc::Halfedge &he) {

  // Note: the missing contribution from faces only contributes to z -
  // axis forces
  // volGrad = vpg->faceNormals[fID] * vpg->faceAreas[fID] /
  // 3;
  std::size_t fID = he.face().getIndex();
  bool interiorHalfedge = he.isInterior();
  gc::Vector3 volGrad{0, 0, 0};
  if (interiorHalfedge) {
    volGrad = vpg.faceNormals[fID] * vpg.faceAreas[fID] / 3;
  }
  return volGrad;
}

gc::VertexData<gc::Vector3> System::computeVertexSchlafliVector() {
  mesh->compress();
  gc::VertexData<gc::Vector3> vector(*mesh, {0, 0, 0});
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};
    double Hi = vpg->vertexMeanCurvatures[i] / vpg->vertexDualAreas[i];
    double H0i = H0[i];
    for (gc::Halfedge he : v.outgoingHalfedges()) {
      std::size_t i_vj = he.tipVertex().getIndex();
      double Hj = vpg->vertexMeanCurvatures[i_vj] / vpg->vertexDualAreas[i_vj];
      double H0j = H0[i_vj];
      gc::Vector3 vec1;
      gc::Vector3 vec2;
      std::tie(vec1, vec2) = computeHalfedgeSchlafliVector(*vpg, he);
      vector[v] += (Hi - H0i) * vec1 + (Hj - H0j) * vec2;
    }
  }
  return vector;
}

gc::VertexData<gc::Vector3> System::computeVertexGaussianCurvatureVector() {
  return halfedgeVectorToVertexVector(*mesh, *vpg,
                                      computeHalfedgeGaussianCurvatureVector);
}

gc::VertexData<gc::Vector3> System::computeVertexMeanCurvatureVector() {
  return halfedgeVectorToVertexVector(*mesh, *vpg,
                                      computeHalfedgeMeanCurvatureVector);
}

gc::VertexData<gc::Vector3> System::computeVertexVolumeVariationVector() {
  return halfedgeVectorToVertexVector(*mesh, *vpg,
                                      computeHalfedgeVolumeVariationVector);
}

gcs::VertexData<gc::Vector3> System::halfedgeVectorToVertexVector(
    gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &vpg,
    std::function<gc::Vector3(gcs::VertexPositionGeometry &vpg, gc::Halfedge &)>
        computeHalfedgeVariationalVector) {
  mesh.compress();
  gc::VertexData<gc::Vector3> vector(mesh, {0, 0, 0});
  for (std::size_t i = 0; i < mesh.nVertices(); ++i) {
    gc::Vertex v{mesh.vertex(i)};
    for (gc::Halfedge he : v.outgoingHalfedges()) {
      vector[v] += computeHalfedgeVariationalVector(vpg, he);
    }
  }
  return vector;
}

void System::computeMechanicalForces() {
  assert(mesh->isCompressed());
  // if(!mesh->isCompressed()){
  //   mem3dg_runtime_error("Mesh must be compressed to compute forces!");
  // }

  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    computeMechanicalForces(i);
  }

  // measure smoothness
  // if (meshProcessor.meshMutator.isSplitEdge ||
  //     meshProcessor.meshMutator.isCollapseEdge) {
  //   isSmooth = !hasOutlier(forces.bendingForce.raw());
  // }
}

void System::computeMechanicalForces(gcs::Vertex &v) {
  size_t i = v.getIndex();
  computeMechanicalForces(i);
}

void System::computeMechanicalForces(size_t i) {
  gc::Vertex v{mesh->vertex(i)};
  gc::Vector3 bendingForceVec{0, 0, 0};
  gc::Vector3 bendingForceVec_areaGrad{0, 0, 0};
  gc::Vector3 bendingForceVec_gaussVec{0, 0, 0};
  gc::Vector3 bendingForceVec_schlafliVec{0, 0, 0};

  gc::Vector3 deviatoricForceVec{0, 0, 0};
  gc::Vector3 deviatoricForceVec_mean{0, 0, 0};
  gc::Vector3 deviatoricForceVec_gauss{0, 0, 0};

  gc::Vector3 capillaryForceVec{0, 0, 0};
  gc::Vector3 osmoticForceVec{0, 0, 0};
  gc::Vector3 lineCapForceVec{0, 0, 0};
  gc::Vector3 adsorptionForceVec{0, 0, 0};
  gc::Vector3 aggregationForceVec{0, 0, 0};
  double Hi = vpg->vertexMeanCurvatures[i] / vpg->vertexDualAreas[i];
  double H0i = H0[i];
  double Kbi = Kb[i];
  double Kdi = Kd[i];
  double proteinDensityi = proteinDensity[i];
  bool boundaryVertex = v.isBoundary();

  for (gc::Halfedge he : v.outgoingHalfedges()) {
    std::size_t fID = he.face().getIndex();

    // Initialize local variables for computation
    std::size_t i_vj = he.tipVertex().getIndex();

    gc::Vector3 dphi_ijk{he.isInterior() ? proteinDensityGradient[fID]
                                         : gc::Vector3{0, 0, 0}};
    double Hj = vpg->vertexMeanCurvatures[i_vj] / vpg->vertexDualAreas[i_vj];
    double H0j = H0[i_vj];
    double Kbj = Kb[i_vj];
    double Kdj = Kd[i_vj];
    double proteinDensityj = proteinDensity[i_vj];
    bool interiorHalfedge = he.isInterior();
    bool boundaryEdge = he.edge().isBoundary();
    bool boundaryNeighborVertex = he.next().vertex().isBoundary();

    gc::Vector3 areaGrad = 2 * computeHalfedgeMeanCurvatureVector(*vpg, he);
    gc::Vector3 gaussVec = computeHalfedgeGaussianCurvatureVector(*vpg, he);
    gc::Vector3 schlafliVec1;
    gc::Vector3 schlafliVec2;
    // std::tie(schlafliVec1, schlafliVec2) =
    //     computeHalfedgeSchlafliVector(*vpg, he);
    schlafliVec1 =
        vpg->edgeLengths[he.edge()] * dihedralAngleGradient(he, he.vertex());
    schlafliVec2 =
        vpg->edgeLengths[he.twin().edge()] *
            dihedralAngleGradient(he.twin(), he.vertex()) +
        vpg->edgeLengths[he.next().edge()] *
            dihedralAngleGradient(he.next(), he.vertex()) +
        vpg->edgeLengths[he.twin().next().next().edge()] *
            dihedralAngleGradient(he.twin().next().next(), he.vertex());
    gc::Vector3 oneSidedAreaGrad{0, 0, 0};
    gc::Vector3 dirichletVec{0, 0, 0};
    if (interiorHalfedge) {
      oneSidedAreaGrad = 0.5 * gc::cross(vpg->faceNormals[fID],
                                         vecFromHalfedge(he.next(), *vpg));
      dirichletVec = computeGradientNorm2Gradient(he, proteinDensity) /
                     vpg->faceAreas[fID];
    }

    // Assemble to forces
    osmoticForceVec +=
        forces.osmoticPressure * computeHalfedgeVolumeVariationVector(*vpg, he);
    capillaryForceVec -= forces.surfaceTension * areaGrad;
    adsorptionForceVec -= (proteinDensityi / 3 + proteinDensityj * 2 / 3) *
                          parameters.adsorption.epsilon * areaGrad;
    aggregationForceVec -= (proteinDensityi * proteinDensityi / 3 +
                            proteinDensityj * proteinDensityj * 2 / 3) *
                           parameters.aggregation.chi * areaGrad;
    lineCapForceVec -=
        parameters.dirichlet.eta *
        (0.125 * dirichletVec - 0.5 * dphi_ijk.norm2() * oneSidedAreaGrad);

    bendingForceVec_schlafliVec -=
        (Kbi * (Hi - H0i) * schlafliVec1 + Kbj * (Hj - H0j) * schlafliVec2);
    bendingForceVec_areaGrad -= (Kbi * (H0i * H0i - Hi * Hi) / 3 +
                                 Kbj * (H0j * H0j - Hj * Hj) * 2 / 3) *
                                areaGrad;
    bendingForceVec_gaussVec -=
        (Kbi * (Hi - H0i) + Kbj * (Hj - H0j)) * gaussVec;

    deviatoricForceVec_mean -=
        (Kdi * Hi + Kdj * Hj) * gaussVec +
        (Kdi * (-Hi * Hi) / 3 + Kdj * (-Hj * Hj) * 2 / 3) * areaGrad +
        (Kdi * Hi * schlafliVec1 + Kdj * Hj * schlafliVec2);

    // bool interiorTwinHalfedge = he.twin().isInterior();
    // if (interiorHalfedge) {
    //   deviatoricForceVec_gauss -=
    //       Kdi * cornerAngleGradient(he.corner(), he.vertex()) +
    //       Kdj * cornerAngleGradient(he.next().corner(), he.vertex());
    // }
    // if (interiorTwinHalfedge) {
    //   deviatoricForceVec_gauss -=
    //       Kdj * cornerAngleGradient(he.twin().corner(), he.vertex());
    // }
    if (boundaryVertex) {
      if (!boundaryEdge)
        deviatoricForceVec_gauss -=
            Kdj * cornerAngleGradient(he.next().corner(), he.vertex()) +
            Kdj * cornerAngleGradient(he.twin().corner(), he.vertex());
    } else {
      if (boundaryNeighborVertex) {
        deviatoricForceVec_gauss -=
            Kdi * cornerAngleGradient(he.corner(), he.vertex());
      } else {
        deviatoricForceVec_gauss -=
            Kdi * cornerAngleGradient(he.corner(), he.vertex()) +
            Kdj * cornerAngleGradient(he.next().corner(), he.vertex()) +
            Kdj * cornerAngleGradient(he.twin().corner(), he.vertex());
      }
    }
  }

  bendingForceVec = bendingForceVec_areaGrad + bendingForceVec_gaussVec +
                    bendingForceVec_schlafliVec;

  // deviatoricForceVec = deviatoricForceVec_gauss;
  // std::cout << "gauss force: " << deviatoricForceVec_gauss << std::ends;
  deviatoricForceVec = deviatoricForceVec_mean + deviatoricForceVec_gauss;
  // deviatoricForceVec = deviatoricForceVec_mean;

  // masking
  bendingForceVec_areaGrad = forces.maskForce(bendingForceVec_areaGrad, i);
  bendingForceVec_gaussVec = forces.maskForce(bendingForceVec_gaussVec, i);
  bendingForceVec_schlafliVec =
      forces.maskForce(bendingForceVec_schlafliVec, i);
  bendingForceVec = forces.maskForce(bendingForceVec, i);

  deviatoricForceVec_mean = forces.maskForce(deviatoricForceVec_mean, i);
  deviatoricForceVec_gauss = forces.maskForce(deviatoricForceVec_gauss, i);
  deviatoricForceVec = forces.maskForce(deviatoricForceVec, i);

  osmoticForceVec = forces.maskForce(osmoticForceVec, i);
  capillaryForceVec = forces.maskForce(capillaryForceVec, i);
  lineCapForceVec = forces.maskForce(lineCapForceVec, i);
  adsorptionForceVec = forces.maskForce(adsorptionForceVec, i);
  aggregationForceVec = forces.maskForce(aggregationForceVec, i);

  // Combine to one
  forces.bendingForceVec_areaGrad[i] = bendingForceVec_areaGrad;
  forces.bendingForceVec_gaussVec[i] = bendingForceVec_gaussVec;
  forces.bendingForceVec_schlafliVec[i] = bendingForceVec_schlafliVec;
  forces.bendingForceVec[i] = bendingForceVec;

  forces.deviatoricForceVec[i] = deviatoricForceVec;
  forces.deviatoricForceVec_mean[i] = deviatoricForceVec_mean;
  forces.deviatoricForceVec_gauss[i] = deviatoricForceVec_gauss;

  forces.capillaryForceVec[i] = capillaryForceVec;
  forces.osmoticForceVec[i] = osmoticForceVec;
  forces.lineCapillaryForceVec[i] = lineCapForceVec;
  forces.adsorptionForceVec[i] = adsorptionForceVec;
  forces.aggregationForceVec[i] = aggregationForceVec;

  // Scalar force by projection to angle-weighted normal
  forces.bendingForce[i] = forces.ontoNormal(bendingForceVec, i);
  forces.deviatoricForce[i] = forces.ontoNormal(deviatoricForceVec, i);
  forces.capillaryForce[i] = forces.ontoNormal(capillaryForceVec, i);
  forces.osmoticForce[i] = forces.ontoNormal(osmoticForceVec, i);
  forces.lineCapillaryForce[i] = forces.ontoNormal(lineCapForceVec, i);
  forces.adsorptionForce[i] = forces.ontoNormal(adsorptionForceVec, i);
  forces.aggregationForce[i] = forces.ontoNormal(aggregationForceVec, i);
}

EigenVectorX3dr System::prescribeExternalForce() {
#define MODE 1
#if MODE == 0 // axial sinusoidal force
  double freq = 5;
  double totalHeight = toMatrix(vpg->inputVertexPositions).col(2).maxCoeff() -
                       toMatrix(vpg->inputVertexPositions).col(2).minCoeff();
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};
    gc::Vector3 direction{vpg->inputVertexPositions[v].x,
                          vpg->inputVertexPositions[v].y, 0};

    double externalPressureMagnitude =
        parameters.external.Kf *
        (1 + sin(freq * 2 * constants::PI / totalHeight *
                 vpg->inputVertexPositions[v].z));
    forces.externalForceVec[i] =
        forces.maskForce(externalPressureMagnitude * vpg->vertexDualArea(v) *
                             direction.normalize(),
                         i);
  }

#elif MODE == 1 // anchor force
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};
    forces.externalForceVec[i] = forces.maskForce(
        parameters.external.Kf *
            ((vpg->vertexGaussianCurvatures[v] < -700 * vpg->vertexDualAreas[v])
                 ? vpg->vertexGaussianCurvatures[v]
                 : 0) *
            vpg->vertexDualAreas[v] * vpg->vertexNormals[v],
        i);
  }

#elif MODE == 2 // anchor force
  double decayTime = 500;
  gcs::HeatMethodDistanceSolver heatSolver(*vpg);
  geodesicDistanceFromPtInd = heatSolver.computeDistance(thePoint);
  double standardDeviation = 0.02;

  // gc::Vector3 anchor{0, 0, 1};
  gc::Vector3 direction{0, 0, 1};
  // direction = anchor - vpg->inputVertexPositions[thePoint.nearestVertex()];
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};
    forces.externalForceVec[i] =
        forces.maskForce(exp(-time / decayTime) * parameters.external.Kf *
                             gaussianDistribution(geodesicDistanceFromPtInd[v],
                                                  standardDeviation) *
                             vpg->vertexDualArea(v) * direction,
                         i);
  }
#endif
  forces.externalForce = forces.ontoNormal(forces.externalForceVec);

  return toMatrix(forces.externalForceVec);
}

void System::computeSelfAvoidanceForce() {
  forces.selfAvoidanceForceVec.fill({0, 0, 0});
  const double d0 = parameters.selfAvoidance.d;
  const double mu = parameters.selfAvoidance.mu;
  const double n = parameters.selfAvoidance.n;
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex vi{mesh->vertex(i)};
    gc::VertexData<bool> neighborList(*mesh, false);
    meshProcessor.meshMutator.markVertices(neighborList, vi, n);
    for (std::size_t j = i + 1; j < mesh->nVertices(); ++j) {
      if (neighborList[j])
        continue;
      gc::Vertex vj{mesh->vertex(j)};
      // double penalty = mu * vpg->vertexDualAreas[vi] * proteinDensity[vi] *
      //                  vpg->vertexDualAreas[vj] * proteinDensity[vj];
      double penalty = mu * proteinDensity[vi] * proteinDensity[vj];
      // double penalty = mu;
      // double penalty = mu * vpg->vertexDualAreas[vi] *
      // vpg->vertexDualAreas[vj];;
      gc::Vector3 r =
          vpg->inputVertexPositions[vj] - vpg->inputVertexPositions[vi];
      double distance = gc::norm(r) - d0;
      gc::Vector3 grad = r.normalize();
      // forces.selfAvoidanceForceVec[i] -=
      //     forces.maskForce(penalty / distance * grad, i);
      // forces.selfAvoidanceForceVec[j] +=
      //     forces.maskForce(penalty / distance * grad, j);
      forces.selfAvoidanceForceVec[i] -=
          forces.maskForce(penalty / distance / distance * grad, i);
      forces.selfAvoidanceForceVec[j] +=
          forces.maskForce(penalty / distance / distance * grad, j);
    }
  }
  forces.selfAvoidanceForce = forces.ontoNormal(forces.selfAvoidanceForceVec);
}

void System::computeChemicalPotentials() {
  gcs::VertexData<double> dH0dphi(*mesh, 0);
  gcs::VertexData<double> dKbdphi(*mesh, 0);
  gcs::VertexData<double> dKddphi(*mesh, 0);
  auto meanCurvDiff = (vpg->vertexMeanCurvatures.raw().array() /
                       vpg->vertexDualAreas.raw().array()) -
                      H0.raw().array();

  if (parameters.bending.relation == "linear") {
    dH0dphi.fill(parameters.bending.H0c);
    dKbdphi.fill(parameters.bending.Kbc);
    dKddphi.fill(parameters.bending.Kdc);
  } else if (parameters.bending.relation == "hill") {
    EigenVectorX1d proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    dH0dphi.raw() =
        (2 * parameters.bending.H0c * proteinDensity.raw().array() /
         ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
            .matrix();
    dKbdphi.raw() =
        (2 * parameters.bending.Kbc * proteinDensity.raw().array() /
         ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
            .matrix();
    dKddphi.raw() =
        (2 * parameters.bending.Kdc * proteinDensity.raw().array() /
         ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
            .matrix();
  }

  forces.bendingPotential.raw() = forces.maskProtein(
      -vpg->vertexDualAreas.raw().array() *
      (meanCurvDiff * meanCurvDiff * dKbdphi.raw().array() -
       2 * Kb.raw().array() * meanCurvDiff * dH0dphi.raw().array()));

  if (parameters.bending.Kd != 0 || parameters.bending.Kdc != 0) {
    forces.deviatoricPotential.raw() =
        -dKddphi.raw().array() *
        (vpg->vertexMeanCurvatures.raw().array().square() /
             vpg->vertexDualAreas.raw().array() -
         vpg->vertexGaussianCurvatures.raw().array());
  }

  if (parameters.adsorption.epsilon != 0)
    forces.adsorptionPotential.raw() = forces.maskProtein(
        -parameters.adsorption.epsilon * vpg->vertexDualAreas.raw().array());

  if (parameters.aggregation.chi != 0)
    forces.aggregationPotential.raw() = forces.maskProtein(
        -2 * parameters.aggregation.chi * proteinDensity.raw().array() *
        vpg->vertexDualAreas.raw().array());

  // if (parameters.adsorption.epsilon != 0)
  //   forces.adsorptionPotential.raw() = forces.maskProtein(
  //       -parameters.adsorption.epsilon * vpg->vertexDualAreas.raw().array() /
  //       vpg->vertexDualAreas.raw().array());

  // if (parameters.aggregation.chi != 0)
  //   forces.aggregationPotential.raw() = forces.maskProtein(
  //       -2 * parameters.aggregation.chi * proteinDensity.raw().array());

  if (parameters.dirichlet.eta != 0)
    forces.diffusionPotential.raw() = forces.maskProtein(
        -parameters.dirichlet.eta * vpg->cotanLaplacian * proteinDensity.raw());

  if (parameters.proteinDistribution.lambdaPhi != 0)
    forces.interiorPenaltyPotential.raw() =
        forces.maskProtein(parameters.proteinDistribution.lambdaPhi *
                           (1 / proteinDensity.raw().array() -
                            1 / (1 - proteinDensity.raw().array())));
  // F.chemicalPotential.raw().array() =
  //     -vpg->vertexDualAreas.raw().array() *
  //     (P.adsorption.epsilon - 2 * Kb.raw().array() * meanCurvDiff *
  //     dH0dphi.raw().array() +
  //      meanCurvDiff * meanCurvDiff * dKbdphi.raw().array());
  // F.chemicalPotential.raw().array() +=
  //     P.proteinDistribution.lambdaPhi * (1 / proteinDensity.raw().array() -
  //                    1 / (1 - proteinDensity.raw().array()));
}

void System::computeDPDForces(double dt) {
  toMatrix(forces.dampingForceVec).setZero();
  toMatrix(forces.stochasticForceVec).setZero();
  // std::default_random_engine random_generator;
  // gcs::EdgeData<double> random_var(mesh);
  double sigma = sqrt(2 * parameters.dpd.gamma * mem3dg::constants::kBoltzmann *
                      parameters.temperature / dt);
  std::normal_distribution<double> normal_dist(0, sigma);

  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    gcs::Vertex v1 = he.vertex();
    gcs::Vertex v2 = he.next().vertex();

    gc::Vector3 dVel12 = velocity[v1] - velocity[v2];
    gc::Vector3 direction =
        (vpg->inputVertexPositions[v1] - vpg->inputVertexPositions[v2])
            .normalize();
    // gc::Vector3 direction =
    //     (vpg->vertexNormals[v1] + vpg->vertexNormals[v2]).normalize();

    gc::Vector3 df =
        parameters.dpd.gamma * (gc::dot(dVel12, direction) * direction);
    forces.dampingForceVec[v1] -= df;
    forces.dampingForceVec[v2] += df;

    if (sigma != 0) {
      double noise = normal_dist(rng);
      forces.stochasticForceVec[v1] += noise * direction;
      forces.stochasticForceVec[v2] -= noise * direction;
    }

    // gc::Vector3 dVel21 = vel[v2] - vel[v1];
    // gc::Vector3 dPos21_n = (pos[v2] - pos[v1]).normalize();

    // std::cout << -gamma * (gc::dot(dVel12, dPos12_n) * dPos12_n)
    //           << " == " << -gamma * (gc::dot(-dVel12, -dPos12_n) *
    //           -dPos12_n)
    //           << " == " << -gamma * (gc::dot(dVel21, dPos21_n) * dPos21_n)
    //           << std::endl;
  }
  forces.dampingForceVec = forces.maskForce(forces.dampingForceVec);
  forces.stochasticForceVec = forces.maskForce(forces.stochasticForceVec);
  // dampingForce_e =
  //     forces.maskForce(forces.addNormal(forces.ontoNormal(dampingForce_e)));
  // stochasticForce_e =
  //     forces.maskForce(forces.addNormal(forces.ontoNormal(stochasticForce_e)));
}

gc::VertexData<gc::Vector3> System::computeDampingForce() {
  return -parameters.damping * velocity;
}

gc::Vector3 System::computeGradientNorm2Gradient(
    const gcs::Halfedge &he, const gcs::VertexData<double> &quantities) {
  if (!he.isInterior()) {
    throw std::runtime_error(
        "computeGradientNormGradient: halfedge is not interior!");
  }

  // quantities
  double qj = quantities[he.next().next().vertex()];
  double qi = quantities[he.vertex()];
  double qk = quantities[he.next().vertex()];

  if (qj == qi && qj == qk) {
    return gc::Vector3({0, 0, 0});
  } else {
    // Edge and normal vector
    gc::Vector3 n = vpg->faceNormals[he.face()];
    gc::Vector3 ej = vecFromHalfedge(he, *vpg);
    gc::Vector3 ei = vecFromHalfedge(he.next(), *vpg);
    gc::Vector3 ek = vecFromHalfedge(he.next().next(), *vpg);

    // exterior angle of triangles (angles formed by e_perp)
    double anglek = gc::angle(ej, ei);
    double anglej = gc::angle(ei, ek);
    double anglei = gc::angle(ek, ej);

    // gradient of edge length wrt he.vertex()
    gc::Vector3 grad_ejnorm = -ej.normalize();
    gc::Vector3 grad_eknorm = ek.normalize();

    // gradient of exterior angle wrt he.vertex()
    gc::Vector3 grad_anglek =
        -cornerAngleGradient(he.next().corner(), he.vertex());
    gc::Vector3 grad_anglej =
        -cornerAngleGradient(he.next().next().corner(), he.vertex());
    gc::Vector3 grad_anglei = -cornerAngleGradient(he.corner(), he.vertex());
    // gc::Vector3 grad_anglek = gc::cross(n, ej).normalize() / gc::norm(ej);
    // gc::Vector3 grad_anglej = gc::cross(n, ek).normalize() / gc::norm(ek);
    // gc::Vector3 grad_anglei = -(grad_anglek + grad_anglej);

    // chain rule
    gc::Vector3 grad_cosanglek = -sin(anglek) * grad_anglek;
    gc::Vector3 grad_cosanglei = -sin(anglei) * grad_anglei;
    gc::Vector3 grad_cosanglej = -sin(anglej) * grad_anglej;

    // g = qj * ej_perp +  qi * ei_perp +  qk * ek_perp
    // gradient of |g|^2
    return 2 * qj * qj * gc::norm(ej) * grad_ejnorm +
           2 * qk * qk * gc::norm(ek) * grad_eknorm +
           2 * qj * qi * gc::norm(ei) *
               (grad_ejnorm * cos(anglek) + gc::norm(ej) * grad_cosanglek) +
           2 * qi * qk * gc::norm(ei) *
               (grad_eknorm * cos(anglej) + gc::norm(ek) * grad_cosanglej) +
           2 * qj * qk *
               (grad_ejnorm * gc::norm(ek) * cos(anglei) +
                gc::norm(ej) * grad_eknorm * cos(anglei) +
                gc::norm(ej) * gc::norm(ek) * grad_cosanglei);
  }
}

double System::computeNorm(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &force) const {
  // Eigen::Matrix<double, Eigen::Dynamic, 1> mask =
  //     outlierMask(force).cast<double>();

  // return rowwiseProduct(mask, force).cwiseAbs().sum() /
  //        rowwiseProduct(mask, vpg->vertexDualAreas.raw()).sum();

  // return force.lpNorm<1>();
  return force.norm();
}

double System::computeNorm(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &&force) const {
  // Eigen::Matrix<double, Eigen::Dynamic, 1> mask =
  //     outlierMask(force).cast<double>();

  // return rowwiseProduct(mask, force).cwiseAbs().sum() /
  //        rowwiseProduct(mask, vpg->vertexDualAreas.raw()).sum();
  // L1 Norm
  // return force.lpNorm<1>();

  return force.norm();
}

void System::computePhysicalForcing() {

  // zero all forces
  forces.mechanicalForceVec.fill({0, 0, 0});

  forces.bendingForceVec.fill({0, 0, 0});
  forces.bendingForceVec_areaGrad.fill({0, 0, 0});
  forces.bendingForceVec_gaussVec.fill({0, 0, 0});
  forces.bendingForceVec_schlafliVec.fill({0, 0, 0});

  forces.deviatoricForceVec.fill({0, 0, 0});

  forces.capillaryForceVec.fill({0, 0, 0});
  forces.osmoticForceVec.fill({0, 0, 0});
  forces.lineCapillaryForceVec.fill({0, 0, 0});
  forces.adsorptionForceVec.fill({0, 0, 0});
  forces.aggregationForceVec.fill({0, 0, 0});
  forces.externalForceVec.fill({0, 0, 0});
  forces.selfAvoidanceForceVec.fill({0, 0, 0});

  forces.dampingForceVec.fill({0, 0, 0});
  forces.stochasticForceVec.fill({0, 0, 0});

  forces.mechanicalForce.raw().setZero();
  forces.bendingForce.raw().setZero();
  forces.deviatoricForce.raw().setZero();
  forces.capillaryForce.raw().setZero();
  forces.lineCapillaryForce.raw().setZero();
  forces.externalForce.raw().setZero();
  forces.adsorptionForce.raw().setZero();
  forces.aggregationForce.raw().setZero();
  forces.osmoticForce.raw().setZero();
  forces.selfAvoidanceForce.raw().setZero();

  forces.chemicalPotential.raw().setZero();

  forces.diffusionPotential.raw().setZero();
  forces.bendingPotential.raw().setZero();
  forces.deviatoricPotential.raw().setZero();
  forces.adsorptionPotential.raw().setZero();
  forces.aggregationPotential.raw().setZero();
  forces.interiorPenaltyPotential.raw().setZero();

  if (parameters.variation.isShapeVariation) {
    computeMechanicalForces();
    if (parameters.external.Kf != 0) {
      prescribeExternalForce();
    }
    if (parameters.selfAvoidance.mu != 0) {
      computeSelfAvoidanceForce();
    }
    forces.mechanicalForceVec =
        forces.osmoticForceVec + forces.capillaryForceVec +
        forces.bendingForceVec + forces.deviatoricForceVec +
        forces.lineCapillaryForceVec + forces.adsorptionForceVec +
        forces.aggregationForceVec + forces.externalForceVec +
        forces.selfAvoidanceForceVec;
    if (parameters.damping != 0)
      forces.mechanicalForceVec += computeDampingForce();
    forces.mechanicalForce = forces.ontoNormal(forces.mechanicalForceVec);
  }

  if (parameters.variation.isProteinVariation) {
    computeChemicalPotentials();
    forces.chemicalPotential =
        forces.adsorptionPotential + forces.aggregationPotential +
        forces.bendingPotential + forces.deviatoricPotential +
        forces.diffusionPotential + forces.interiorPenaltyPotential;
  }

  // compute the mechanical error norm
  mechErrorNorm = parameters.variation.isShapeVariation
                      ? computeNorm(toMatrix(forces.mechanicalForceVec))
                      : 0;

  // compute the chemical error norm
  chemErrorNorm = parameters.variation.isProteinVariation
                      ? computeNorm(forces.chemicalPotential.raw())
                      : 0;
}

void System::computePhysicalForcing(double timeStep) {
  computePhysicalForcing();
  if (parameters.variation.isShapeVariation && parameters.dpd.gamma != 0) {
    computeDPDForces(timeStep);
    forces.mechanicalForceVec +=
        forces.dampingForceVec + forces.stochasticForceVec;
  }

  // if (!f.mesh->hasBoundary()) {
  //   removeTranslation(physicalForceVec);
  //   removeRotation(toMatrix(f.vpg->inputVertexPositions),
  //                  physicalForceVec);
  //   // removeTranslation(DPDPressure);
  //   // removeRotation(toMatrix(f.vpg->inputVertexPositions),
  //   // DPDPressure);
  // }
}

} // namespace solver
} // namespace mem3dg
