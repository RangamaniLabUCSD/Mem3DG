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

gc::Vector3 System::computeMeanCurvatureVector(const gc::Halfedge &he) {
  std::size_t fID = he.face().getIndex();
  std::size_t fID_he_twin = he.twin().face().getIndex();
  bool interiorHalfedge = he.isInterior();
  bool interiorTwinHalfedge = he.twin().isInterior();
  gc::Vector3 areaGrad{0, 0, 0};
  if (interiorHalfedge) {
    areaGrad += 0.25 * gc::cross(vpg->faceNormals[fID],
                                 vecFromHalfedge(he.next(), *vpg));
  }
  if (interiorTwinHalfedge) {
    areaGrad +=
        0.25 * gc::cross(vpg->faceNormals[fID_he_twin],
                         vecFromHalfedge(he.twin().next().next(), *vpg));
  }
  return areaGrad/2;
}

gc::Vector3 System::computeVolumeVariationVector(const gc::Halfedge &he) {
  std::size_t fID = he.face().getIndex();
  bool interiorHalfedge = he.isInterior();
  gc::Vector3 volGrad{0, 0, 0};
  if (interiorHalfedge) {
    volGrad = vpg->faceNormals[fID] * vpg->faceAreas[fID] / 3;
  }
  return volGrad;
}

gc::VertexData<gc::Vector3> System::computeVolumeVariationVector() {
  assert(mesh->isCompressed());
  gc::VertexData<gc::Vector3> vector(*mesh, {0, 0, 0});
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};
    for (gc::Halfedge he : v.outgoingHalfedges()) {
      vector[v] += computeVolumeVariationVector(he);
    }
  }
  return vector;
}

void System::computeVectorForces() {
  assert(mesh->isCompressed());
  // if(!mesh->isCompressed()){
  //   mem3dg_runtime_error("Mesh must be compressed to compute forces!");
  // }

  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};

    gc::Vector3 bendForceVec{0, 0, 0};
    gc::Vector3 bendForceVec_areaGrad{0, 0, 0};
    gc::Vector3 bendForceVec_gaussVec{0, 0, 0};
    gc::Vector3 bendForceVec_schlafliVec{0, 0, 0};

    gc::Vector3 capillaryForceVec{0, 0, 0};
    gc::Vector3 osmoticForceVec{0, 0, 0};
    gc::Vector3 lineCapForceVec{0, 0, 0};
    gc::Vector3 adsorptionForceVec{0, 0, 0};
    double Hi = vpg->vertexMeanCurvatures[i] / vpg->vertexDualAreas[i];
    double H0i = H0[i];
    double Kbi = Kb[i];
    double proteinDensityi = proteinDensity[i];
    bool boundaryVertex = v.isBoundary();

    for (gc::Halfedge he : v.outgoingHalfedges()) {
      std::size_t fID = he.face().getIndex();
      std::size_t fID_he_twin = he.twin().face().getIndex();
      std::size_t heID_twin = he.twin().getIndex();
      std::size_t heID_twin_next = he.twin().next().getIndex();
      std::size_t heID_he_next_next = he.next().next().getIndex();

      // Initialize local variables for computation
      gc::Vertex vj = he.tipVertex();
      std::size_t i_vj = vj.getIndex();

      gc::Vector3 dphi_ijk{he.isInterior() ? proteinDensityGradient[fID]
                                           : gc::Vector3{0, 0, 0}};
      double Hj = vpg->vertexMeanCurvatures[i_vj] / vpg->vertexDualAreas[i_vj];
      double H0j = H0[i_vj];
      double Kbj = Kb[i_vj];
      double proteinDensityj = proteinDensity[i_vj];

      bool boundaryEdge = he.edge().isBoundary();
      bool interiorHalfedge = he.isInterior();
      bool interiorTwinHalfedge = he.twin().isInterior();

      // Note: the missing contribution from faces only contributes to z -
      // axis forces
      // volGrad = vpg->faceNormals[fID] * vpg->faceAreas[fID] /
      // 3;
      gc::Vector3 oneSidedAreaGrad{0, 0, 0};
      gc::Vector3 dirichletVec{0, 0, 0};

      if (interiorHalfedge) {
        oneSidedAreaGrad = 0.5 * gc::cross(vpg->faceNormals[fID],
                                           vecFromHalfedge(he.next(), *vpg));
        dirichletVec = computeGradientNorm2Gradient(he, proteinDensity) /
                       vpg->faceAreas[fID];
      }

      gc::Vector3 areaGrad = 2 * computeMeanCurvatureVector(he);
      gc::Vector3 gaussVec{0, 0, 0};
      gc::Vector3 schlafliVec1{0, 0, 0};
      gc::Vector3 schlafliVec2{0, 0, 0};

      if (!boundaryEdge) {
        // gc::Vector3 eji{} = -vecFromHalfedge(he, *vpg);
        gaussVec = 0.5 * vpg->edgeDihedralAngles[he.edge()] *
                   (-vecFromHalfedge(he, *vpg)).unit();
        schlafliVec1 = vpg->halfedgeCotanWeights[heID_he_next_next] *
                           vpg->faceNormals[fID] +
                       vpg->halfedgeCotanWeights[heID_twin_next] *
                           vpg->faceNormals[fID_he_twin];
      }

      if (boundaryVertex && boundaryEdge) {
        schlafliVec2 = interiorHalfedge
                           ? (-(vpg->halfedgeCotanWeights[he] +
                                vpg->halfedgeCotanWeights[heID_he_next_next]) *
                              vpg->faceNormals[fID])
                           : (-(vpg->halfedgeCotanWeights[heID_twin] +
                                vpg->halfedgeCotanWeights[heID_twin_next]) *
                              vpg->faceNormals[fID_he_twin]);
      } else if (!boundaryVertex && vj.isBoundary()) {
        schlafliVec2 = vpg->halfedgeCotanWeights[heID_he_next_next] *
                           vpg->faceNormals[fID] +
                       vpg->halfedgeCotanWeights[heID_twin_next] *
                           vpg->faceNormals[fID_he_twin];

        if (!he.next().edge().isBoundary())
          schlafliVec2 -= (vpg->halfedgeCotanWeights[he] +
                           vpg->halfedgeCotanWeights[heID_he_next_next]) *
                          vpg->faceNormals[fID];

        if (!he.twin().next().next().edge().isBoundary())
          schlafliVec2 -= (vpg->halfedgeCotanWeights[heID_twin] +
                           vpg->halfedgeCotanWeights[heID_twin_next]) *
                          vpg->faceNormals[fID_he_twin];
      } else {
        schlafliVec2 = -(vpg->halfedgeCotanWeights[he] * vpg->faceNormals[fID] +
                         vpg->halfedgeCotanWeights[heID_twin] *
                             vpg->faceNormals[fID_he_twin]);
      }

      // Assemble to forces
      osmoticForceVec +=
          forces.osmoticPressure * computeVolumeVariationVector(he);
      capillaryForceVec -= forces.surfaceTension * areaGrad;
      adsorptionForceVec -= (proteinDensityi / 3 + proteinDensityj * 2 / 3) *
                            parameters.adsorption.epsilon * areaGrad;
      lineCapForceVec -=
          parameters.dirichlet.eta *
          (0.125 * dirichletVec - 0.5 * dphi_ijk.norm2() * oneSidedAreaGrad);

      bendForceVec_schlafliVec -=
          (Kbi * (Hi - H0i) * schlafliVec1 + Kbj * (Hj - H0j) * schlafliVec2);
      bendForceVec_areaGrad -= (Kbi * (H0i * H0i - Hi * Hi) / 3 +
                                Kbj * (H0j * H0j - Hj * Hj) * 2 / 3) *
                               areaGrad;
      bendForceVec_gaussVec -= (Kbi * (Hi - H0i) + Kbj * (Hj - H0j)) * gaussVec;
      bendForceVec -=
          (Kbi * (Hi - H0i) + Kbj * (Hj - H0j)) * gaussVec +
          (Kbi * (H0i * H0i - Hi * Hi) / 3 +
           Kbj * (H0j * H0j - Hj * Hj) * 2 / 3) *
              areaGrad +
          (Kbi * (Hi - H0i) * schlafliVec1 + Kbj * (Hj - H0j) * schlafliVec2);

      // Compare principal curvature vs truncated curvature vector

      // This is outside the loop
      // gc::Vector3 truncatedCurv{0, 0, 0};
      // gc::Vector3 curvature{0, 0, 0};

      // bool isinterface = false;
      // for (gcs::Face f : v.adjacentFaces()) {
      //   if (dH0[f].norm() > 1) {
      //     isinterface = true;
      //   }
      // }

      // inside of loop
      //   if (isinterface) {
      //     truncatedCurv +=
      //         oneSidedAreaGrad - gc::dot(oneSidedAreaGrad, dH0[fID])
      //         *
      //                                dH0[fID] /
      //                                dH0[fID].norm2();
      //     curvature += oneSidedAreaGrad;
      //   }
      // }

      // if (isinterface) {
      //   std::cout << "curvature: " << curvature.norm()
      //             << " and ontoNormal: " << F.ontoNormal(curvature, v)
      //             << std::endl;
      //   std::cout << "capillaryForceVec / surfacetension: "
      //             << (capillaryForceVec / surfaceTension).norm()
      //             << " and ontoNormal: "
      //             << F.ontoNormal((capillaryForceVec / surfaceTension), v)
      //             << std::endl;
      //   std::cout << "truncated curvature: " << truncatedCurv.norm()
      //             << " and ontoNormal: " << F.ontoNormal(truncatedCurv, v)
      //             << std::endl;
      //   vpg->requireVertexMaxPrincipalCurvatures();
      //   vpg->requireVertexMinPrincipalCurvatures();
      //   std::cout << "principal curvature: "
      //             << vpg->vertexMaxPrincipalCurvatures[v] *
      //                    vpg->vertexDualAreas[v]
      //             << " and: "
      //             << vpg->vertexMinPrincipalCurvatures[v] *
      //                    vpg->vertexDualAreas[v]
      //             << std::endl;
    }

    // masking
    osmoticForceVec = forces.maskForce(osmoticForceVec, i);
    capillaryForceVec = forces.maskForce(capillaryForceVec, i);
    bendForceVec = forces.maskForce(bendForceVec, i);
    lineCapForceVec = forces.maskForce(lineCapForceVec, i);
    adsorptionForceVec = forces.maskForce(adsorptionForceVec, i);
    bendForceVec = forces.maskForce(bendForceVec, i);
    bendForceVec_areaGrad = forces.maskForce(bendForceVec_areaGrad, i);
    bendForceVec_gaussVec = forces.maskForce(bendForceVec_gaussVec, i);
    bendForceVec_schlafliVec = forces.maskForce(bendForceVec_schlafliVec, i);

    // Combine to one
    forces.osmoticForceVec[i] = osmoticForceVec;
    forces.capillaryForceVec[i] = capillaryForceVec;
    forces.bendingForceVec[i] = bendForceVec;
    forces.lineCapillaryForceVec[i] = lineCapForceVec;
    forces.adsorptionForceVec[i] = adsorptionForceVec;
    forces.bendingForceVec[i] = bendForceVec;
    forces.bendingForceVec_areaGrad[i] = bendForceVec_areaGrad;
    forces.bendingForceVec_gaussVec[i] = bendForceVec_gaussVec;
    forces.bendingForceVec_schlafliVec[i] = bendForceVec_schlafliVec;

    // Scalar force by projection to angle-weighted normal
    forces.bendingForce[i] = forces.ontoNormal(bendForceVec, i);
    forces.capillaryForce[i] = forces.ontoNormal(capillaryForceVec, i);
    forces.osmoticForce[i] = forces.ontoNormal(osmoticForceVec, i);
    forces.lineCapillaryForce[i] = forces.ontoNormal(lineCapForceVec, i);
  }

  // measure smoothness
  if (meshProcessor.meshMutator.isSplitEdge ||
      meshProcessor.meshMutator.isCollapseEdge) {
    isSmooth = !hasOutlier(forces.bendingForce.raw());
  }
}

EigenVectorX1d System::computeBendingForce() {
  mem3dg_runtime_error("Out of data implementation, shouldn't be called!");
  // A. non-optimized version
  // if (O.isLocalCurvature) {
  //   // Split calculation for two domain
  //   bendingPressure.raw().setZero();
  //   auto subdomain = [&](double H0_temp) {
  //     EigenVectorX1d lap_H = vpg->vertexLumpedMassMatrix.cwiseInverse() * L
  //     * (H.raw().array() - H0_temp).matrix(); EigenVectorX1d scalerTerms =
  //         rowwiseProduct(H.raw(), H.raw()) + H.raw() * H0_temp - K.raw();
  //     EigenVectorX1d productTerms =
  //         2.0 *
  //         rowwiseProduct(scalerTerms, (H.raw().array() -
  //         H0_temp).matrix());
  //     bendingPressure.raw().array() +=
  //         (H0.raw().array() == H0_temp).cast<double>().array() *
  //         (-P.Kb * (productTerms + lap_H)).array();
  //   };
  //   subdomain(P.H0);
  //   subdomain(0);
  // } else {

  EigenVectorX1d ptwiseH = vpg->vertexMeanCurvatures.raw().array() /
                           vpg->vertexDualAreas.raw().array();

  // calculate the Laplacian of mean curvature H
  EigenVectorX1d lap_H =
      -(vpg->cotanLaplacian * Kb.raw().cwiseProduct(ptwiseH - H0.raw()))
           .array() /
      vpg->vertexDualAreas.raw().array();

  // initialize and calculate intermediary result scalerTerms
  EigenVectorX1d scalerTerms = ptwiseH.cwiseProduct(ptwiseH) +
                               ptwiseH.cwiseProduct(H0.raw()) -
                               (vpg->vertexGaussianCurvatures.raw().array() /
                                vpg->vertexDualAreas.raw().array())
                                   .matrix();
  // scalerTerms = scalerTerms.array().max(0);

  // initialize and calculate intermediary result productTerms
  EigenVectorX1d productTerms =
      -2.0 *
      (Kb.raw().array() * (ptwiseH - H0.raw()).array() * scalerTerms.array())
          .matrix();

  // calculate bendingForce
  forces.bendingForce.raw() =
      vpg->vertexLumpedMassMatrix * (productTerms + lap_H);
  // }

  isSmooth = !hasOutlier(forces.bendingForce.raw());

  return forces.bendingForce.raw();

  // /// B. optimized version
  // // calculate the Laplacian of mean curvature H
  // EigenVectorX1d lap_H_integrated = L * (H - H0);

  // // initialize and calculate intermediary result scalarTerms_integrated
  // EigenVectorX1d H_integrated = M * H;
  // EigenVectorX1d scalarTerms_integrated =
  //     M * rowwiseProduct(vpg->vertexLumpedMassMatrix.cwiseInverse() *
  //     H_integrated, vpg->vertexLumpedMassMatrix.cwiseInverse() *
  //     H_integrated) + rowwiseProduct(H_integrated, H0) -
  //     vpg->vertexGaussianCurvatures.raw();
  // EigenVectorX1d zeroMatrix;
  // zeroMatrix.resize(n_vertices, 1);
  // zeroMatrix.setZero();
  // scalarTerms_integrated =
  //     scalarTerms_integrated.array().max(zeroMatrix.array());

  // // initialize and calculate intermediary result productTerms_integrated
  // EigenVectorX1d productTerms_integrated;
  // productTerms_integrated.resize(n_vertices, 1);
  // productTerms_integrated =
  //     2.0 * rowwiseProduct(scalarTerms_integrated, H - H0);

  // bendingPressure_e =
  //     -2.0 * P.Kb *
  //     rowwiseScalarProduct(vpg->vertexLumpedMassMatrix.cwiseInverse() *
  //     (productTerms_integrated + lap_H_integrated),
  //                    vertexAngleNormal_e);
}

EigenVectorX1d System::computeCapillaryForce() {
  mem3dg_runtime_error("computeCapillaryForce: out of data implementation, "
                       "shouldn't be called!");
  /// Geometric implementation
  forces.surfaceTension = parameters.tension.Ksg *
                              (surfaceArea - parameters.tension.At) /
                              parameters.tension.At +
                          parameters.tension.lambdaSG;
  forces.capillaryForce.raw() =
      -forces.surfaceTension * 2 * vpg->vertexMeanCurvatures.raw();

  return forces.capillaryForce.raw();

  // /// Nongeometric implementationx
  // for (gcs::Vertex v : mesh->vertices()) {
  //   gc::Vector3 globalForce{0.0, 0.0, 0.0};
  //   for (gcs::Halfedge he : v.outgoingHalfedges()) {
  //     gc::Vector3 base_vec = vecFromHalfedge(he.next(), vpg);
  //     gc::Vector3 localAreaGradient =
  //         -gc::cross(base_vec, vpg->faceNormals[he.face()]);
  //     assert((gc::dot(localAreaGradient, vecFromHalfedge(he, vpg))) < 0);
  //     if (P.Ksg != 0) {
  //       capillaryPressure[v] += -P.Ksg * localAreaGradient *
  //                               (surfaceArea - refSurfaceArea) /
  //                               refSurfaceArea;
  //     }
  //   }
  //   capillaryPressure[v] /= vpg->vertexDualAreas[v];
  // }
}

EigenVectorX1d System::computeOsmoticForce() {
  mem3dg_runtime_error("computeOsmoticForce: out of data implementation, "
                       "shouldn't be called!");
  forces.osmoticForce.raw().setConstant(forces.osmoticPressure);
  forces.osmoticForce.raw().array() *= vpg->vertexDualAreas.raw().array();

  return forces.osmoticForce.raw();

  // /// Nongeometric implementation
  // for (gcs::Vertex v : mesh->vertices()) {
  //   for (gcs::Halfedge he : v.outgoingHalfedges()) {
  //     gc::Vector3 p1 = vpg->inputVertexPositions[he.next().vertex()];
  //     gc::Vector3 p2 =
  //     vpg->inputVertexPositions[he.next().next().vertex()]; gc::Vector3
  //     dVdx = 0.5 * gc::cross(p1, p2) / 3.0; insidePressure[v] +=
  //         -P.Kv * (volume - refVolume * P.Vt) / (refVolume * P.Vt) * dVdx;
  //   }
  // }
}

EigenVectorX1d System::computeLineCapillaryForce() {
  if (false) {
    mem3dg_runtime_error(
        "computeLineCapillaryForce: out of data implementation, "
        "shouldn't be called!");
    // zeros out the nonpositive normal curvature to compensate the fact that
    // d0 is ill-defined in low resolution auto normalCurvature =
    // vpg->edgeDihedralAngles.raw(); F.lineCapillaryForce.raw() =
    //     -D * vpg->hodge1Inverse *
    //     ((vpg->hodge1 *
    //       (F.lineTension.raw().array() / vpg->edgeLengths.raw().array())
    //           .matrix())
    //          .array() *
    //      normalCurvature.array().max(0))
    //         .matrix();
  }
  return forces.lineCapillaryForce.raw();
}

EigenVectorX1d System::computeExternalForce() {
#define MODE 0
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
    forces.externalForceVec[i] = externalPressureMagnitude *
                                 vpg->vertexDualArea(v) * direction.normalize();
  }

#elif MODE == 1 // anchor force
  double concentration = 10;
  gcs::HeatMethodDistanceSolver heatSolver(*vpg);
  geodesicDistanceFromPtInd = heatSolver.computeDistance(thePoint);
  double standardDeviation =
      geodesicDistanceFromPtInd.raw().maxCoeff() / concentration;

  gc::Vector3 anchor{0, 0, 1};
  gc::Vector3 direction;
  direction = anchor - vpg->inputVertexPositions[thePoint.nearestVertex()];
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};
    forces.externalForceVec[i] =
        parameters.external.Kf *
        gaussianDistribution(geodesicDistanceFromPtInd[v], standardDeviation) *
        vpg->vertexDualArea(v) * direction;
  }
#endif
  forces.externalForce = forces.ontoNormal(forces.externalForceVec);
  return forces.externalForce.raw();
}

EigenVectorX1d System::computeChemicalPotential() {
  gcs::VertexData<double> dH0dphi(*mesh, 0);
  gcs::VertexData<double> dKbdphi(*mesh, 0);
  auto meanCurvDiff = (vpg->vertexMeanCurvatures.raw().array() /
                       vpg->vertexDualAreas.raw().array()) -
                      H0.raw().array();

  if (parameters.bending.relation == "linear") {
    dH0dphi.fill(parameters.bending.H0c);
    dKbdphi.fill(parameters.bending.Kbc);
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
  }

  forces.adsorptionPotential.raw() = forces.maskProtein(
      -parameters.adsorption.epsilon * vpg->vertexDualAreas.raw().array());
  forces.bendingPotential.raw() = forces.maskProtein(
      -vpg->vertexDualAreas.raw().array() *
      (meanCurvDiff * meanCurvDiff * dKbdphi.raw().array() -
       2 * Kb.raw().array() * meanCurvDiff * dH0dphi.raw().array()));
  forces.diffusionPotential.raw() = forces.maskProtein(
      -parameters.dirichlet.eta * vpg->cotanLaplacian * proteinDensity.raw());
  forces.interiorPenaltyPotential.raw() =
      forces.maskProtein(parameters.proteinDistribution.lambdaPhi *
                         (1 / proteinDensity.raw().array() -
                          1 / (1 - proteinDensity.raw().array())));
  forces.chemicalPotential.raw() =
      forces.adsorptionPotential.raw() + forces.bendingPotential.raw() +
      forces.diffusionPotential.raw() + forces.interiorPenaltyPotential.raw();

  // F.chemicalPotential.raw().array() =
  //     -vpg->vertexDualAreas.raw().array() *
  //     (P.adsorption.epsilon - 2 * Kb.raw().array() * meanCurvDiff *
  //     dH0dphi.raw().array() +
  //      meanCurvDiff * meanCurvDiff * dKbdphi.raw().array());
  // F.chemicalPotential.raw().array() +=
  //     P.proteinDistribution.lambdaPhi * (1 / proteinDensity.raw().array() -
  //                    1 / (1 - proteinDensity.raw().array()));

  return forces.chemicalPotential.raw();
}

std::tuple<EigenVectorX3dr, EigenVectorX3dr>
System::computeDPDForces(double dt) {

  auto dampingForce_e = EigenMap<double, 3>(forces.dampingForce);
  auto stochasticForce_e = EigenMap<double, 3>(forces.stochasticForce);

  // Reset forces to zero
  dampingForce_e.setZero();
  stochasticForce_e.setZero();

  // alias positions
  const auto &pos = vpg->inputVertexPositions;

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
    gc::Vector3 dPos12_n = (pos[v1] - pos[v2]).normalize();

    if (parameters.dpd.gamma != 0) {
      gc::Vector3 df =
          parameters.dpd.gamma * (gc::dot(dVel12, dPos12_n) * dPos12_n);
      forces.dampingForce[v1] -= df;
      forces.dampingForce[v2] += df;
    }

    if (sigma != 0) {
      double noise = normal_dist(rng);
      forces.stochasticForce[v1] += noise * dPos12_n;
      forces.stochasticForce[v2] -= noise * dPos12_n;
    }

    // gc::Vector3 dVel21 = vel[v2] - vel[v1];
    // gc::Vector3 dPos21_n = (pos[v2] - pos[v1]).normalize();

    // std::cout << -gamma * (gc::dot(dVel12, dPos12_n) * dPos12_n)
    //           << " == " << -gamma * (gc::dot(-dVel12, -dPos12_n) *
    //           -dPos12_n)
    //           << " == " << -gamma * (gc::dot(dVel21, dPos21_n) * dPos21_n)
    //           << std::endl;
  }

  return std::tie(dampingForce_e, stochasticForce_e);
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
    gc::Vector3 grad_anglek = gc::cross(n, ej).normalize() / gc::norm(ej);
    gc::Vector3 grad_anglej = gc::cross(n, ek).normalize() / gc::norm(ek);
    gc::Vector3 grad_anglei = -(grad_anglek + grad_anglej);

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

void System::computePhysicalForces() {

  // zero all forces
  forces.bendingForceVec.fill({0, 0, 0});
  forces.bendingForceVec_areaGrad.fill({0, 0, 0});
  forces.bendingForceVec_gaussVec.fill({0, 0, 0});
  forces.bendingForceVec_schlafliVec.fill({0, 0, 0});

  forces.capillaryForceVec.fill({0, 0, 0});
  forces.osmoticForceVec.fill({0, 0, 0});
  forces.lineCapillaryForceVec.fill({0, 0, 0});
  forces.adsorptionForceVec.fill({0, 0, 0});
  forces.externalForceVec.fill({0, 0, 0});

  forces.bendingForce.raw().setZero();
  forces.capillaryForce.raw().setZero();
  forces.lineCapillaryForce.raw().setZero();
  forces.externalForce.raw().setZero();
  forces.osmoticForce.raw().setZero();

  forces.chemicalPotential.raw().setZero();
  forces.diffusionPotential.raw().setZero();
  forces.bendingPotential.raw().setZero();
  forces.adsorptionPotential.raw().setZero();
  forces.interiorPenaltyPotential.raw().setZero();

  if (parameters.variation.isShapeVariation) {
    computeVectorForces();
    if (parameters.external.Kf != 0) {
      computeExternalForce();
    }
    forces.mechanicalForceVec =
        forces.osmoticForceVec + forces.capillaryForceVec +
        forces.bendingForceVec + forces.lineCapillaryForceVec +
        forces.adsorptionForceVec + forces.externalForceVec;
    forces.mechanicalForce = forces.ontoNormal(forces.mechanicalForceVec);
  }

  if (parameters.variation.isProteinVariation) {
    computeChemicalPotential();
  }

  // computeBendingForce();
  // if (P.Kv != 0) {
  //   computeOsmoticForce();
  // }
  // if (P.Ksg != 0) {
  //   computeCapillaryForce();
  // }
  // if (P.dirichlet.eta != 0) {
  //   computeLineCapillaryForce();
  // }

  // compute the mechanical error norm
  mechErrorNorm = parameters.variation.isShapeVariation
                      ? computeNorm(toMatrix(forces.mechanicalForceVec))
                      : 0;

  // compute the chemical error norm
  chemErrorNorm = parameters.variation.isProteinVariation
                      ? computeNorm(forces.chemicalPotential.raw())
                      : 0;
}

} // namespace solver
} // namespace mem3dg
