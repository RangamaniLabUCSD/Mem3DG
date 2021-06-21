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

void System::computeVectorForces() {
  assert(mesh->isCompressed());
  // if(!mesh->isCompressed()){
  //   mem3dg_runtime_error("Mesh must be compressed to compute forces!");
  // }

  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};

    gc::Vector3 bendForceVec{0, 0, 0};
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
      gc::Vector3 volGrad{0, 0, 0};
      gc::Vector3 oneSidedAreaGrad{0, 0, 0};
      gc::Vector3 dirichletVec{0, 0, 0};
      gc::Vector3 areaGrad{0, 0, 0};

      if (interiorHalfedge) {
        volGrad = vpg->faceNormals[fID] * vpg->faceAreas[fID] / 3;
        oneSidedAreaGrad = 0.5 * gc::cross(vpg->faceNormals[fID],
                                           vecFromHalfedge(he.next(), *vpg));
        dirichletVec = computeGradientNorm2Gradient(he, proteinDensity) /
                       vpg->faceAreas[fID];

        areaGrad = 0.25 * gc::cross(vpg->faceNormals[fID],
                                    vecFromHalfedge(he.next(), *vpg));
        if (interiorTwinHalfedge)
          areaGrad +=
              0.25 * gc::cross(vpg->faceNormals[fID_he_twin],
                               vecFromHalfedge(he.twin().next().next(), *vpg));
      }

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
          schlafliVec2 += -(vpg->halfedgeCotanWeights[he] +
                            vpg->halfedgeCotanWeights[heID_he_next_next]) *
                          vpg->faceNormals[fID];

        if (!he.twin().next().next().edge().isBoundary())
          schlafliVec2 += -(vpg->halfedgeCotanWeights[heID_twin] +
                            vpg->halfedgeCotanWeights[heID_twin_next]) *
                          vpg->faceNormals[fID_he_twin];
      } else {
        schlafliVec2 = -(vpg->halfedgeCotanWeights[he] * vpg->faceNormals[fID] +
                         vpg->halfedgeCotanWeights[heID_twin] *
                             vpg->faceNormals[fID_he_twin]);
      }

      // Assemble to forces
      osmoticForceVec += F.osmoticPressure * volGrad;
      capillaryForceVec -= F.surfaceTension * areaGrad;
      adsorptionForceVec -= (proteinDensityi / 3 + proteinDensityj * 2 / 3) *
                            P.epsilon * areaGrad;
      lineCapForceVec -= P.eta * (0.125 * dirichletVec -
                                  0.5 * dphi_ijk.norm2() * oneSidedAreaGrad);
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
    osmoticForceVec = F.maskForce(osmoticForceVec, i);
    capillaryForceVec = F.maskForce(capillaryForceVec, i);
    bendForceVec = F.maskForce(bendForceVec, i);
    lineCapForceVec = F.maskForce(lineCapForceVec, i);
    adsorptionForceVec = F.maskForce(adsorptionForceVec, i);

    // Combine to one
    F.osmoticForceVec[i] = osmoticForceVec;
    F.capillaryForceVec[i] = capillaryForceVec;
    F.bendingForceVec[i] = bendForceVec;
    F.lineCapillaryForceVec[i] = lineCapForceVec;
    F.adsorptionForceVec[i] = adsorptionForceVec;

    // Scalar force by projection to angle-weighted normal
    F.bendingForce[i] = F.ontoNormal(bendForceVec, i);
    F.capillaryForce[i] = F.ontoNormal(capillaryForceVec, i);
    F.osmoticForce[i] = F.ontoNormal(osmoticForceVec, i);
    F.lineCapillaryForce[i] = F.ontoNormal(lineCapForceVec, i);
  }

  // measure smoothness
  if (O.isSplitEdge || O.isCollapseEdge) {
    isSmooth = !hasOutlier(F.bendingForce.raw());
  }
}

EigenVectorX1d System::computeBendingForce() {
  mem3dg_runtime_error("Out of data implementation, shouldn't be called!");
  // A. non-optimized version
  // if (O.isLocalCurvature) {
  //   // Split calculation for two domain
  //   bendingPressure.raw().setZero();
  //   auto subdomain = [&](double H0_temp) {
  //     EigenVectorX1d lap_H = vpg->vertexLumpedMassMatrix.cwiseInverse() * L *
  //     (H.raw().array() - H0_temp).matrix(); EigenVectorX1d scalerTerms =
  //         rowwiseProduct(H.raw(), H.raw()) + H.raw() * H0_temp - K.raw();
  //     EigenVectorX1d productTerms =
  //         2.0 *
  //         rowwiseProduct(scalerTerms, (H.raw().array() - H0_temp).matrix());
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
  F.bendingForce.raw() = vpg->vertexLumpedMassMatrix * (productTerms + lap_H);
  // }

  isSmooth = !hasOutlier(F.bendingForce.raw());

  return F.bendingForce.raw();

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
  throw std::runtime_error("computeCapillaryForce: out of data implementation, "
                           "shouldn't be called!");
  /// Geometric implementation
  F.surfaceTension =
      P.Ksg * (surfaceArea - refSurfaceArea) / refSurfaceArea + P.lambdaSG;
  F.capillaryForce.raw() =
      -F.surfaceTension * 2 * vpg->vertexMeanCurvatures.raw();

  return F.capillaryForce.raw();

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
  throw std::runtime_error("computeOsmoticForce: out of data implementation, "
                           "shouldn't be called!");
  F.osmoticForce.raw().setConstant(F.osmoticPressure);
  F.osmoticForce.raw().array() *= vpg->vertexDualAreas.raw().array();

  return F.osmoticForce.raw();

  // /// Nongeometric implementation
  // for (gcs::Vertex v : mesh->vertices()) {
  //   for (gcs::Halfedge he : v.outgoingHalfedges()) {
  //     gc::Vector3 p1 = vpg->inputVertexPositions[he.next().vertex()];
  //     gc::Vector3 p2 = vpg->inputVertexPositions[he.next().next().vertex()];
  //     gc::Vector3 dVdx = 0.5 * gc::cross(p1, p2) / 3.0;
  //     insidePressure[v] +=
  //         -P.Kv * (volume - refVolume * P.Vt) / (refVolume * P.Vt) * dVdx;
  //   }
  // }
}

EigenVectorX1d System::computeLineCapillaryForce() {
  if (false) {
    throw std::runtime_error(
        "computeLineCapillaryForce: out of data implementation, "
        "shouldn't be called!");
    // zeros out the nonpositive normal curvature to compensate the fact that d0
    // is ill-defined in low resolution
    // auto normalCurvature = vpg->edgeDihedralAngles.raw();
    // F.lineCapillaryForce.raw() =
    //     -D * vpg->hodge1Inverse *
    //     ((vpg->hodge1 *
    //       (F.lineTension.raw().array() / vpg->edgeLengths.raw().array())
    //           .matrix())
    //          .array() *
    //      normalCurvature.array().max(0))
    //         .matrix();
  }
  return F.lineCapillaryForce.raw();
}

EigenVectorX1d System::computeExternalForce() {
  EigenVectorX1d externalPressureMagnitude;

  // a. FIND OUT THE CURRENT EXTERNAL PRESSURE MAGNITUDE BASED ON CURRENT
  // GEOMETRY

  // auto &dist_e = heatMethodDistance(vpg, mesh->vertex(P.ptInd)).raw();
  // double stdDev = dist_e.maxCoeff() / P.conc;
  // externalPressureMagnitude =
  //    P.Kf / (stdDev * pow(M_PI * 2, 0.5)) *
  //    (-dist_e.array() * dist_e.array() / (2 * stdDev * stdDev)).exp();

  // b. APPLY EXTERNAL PRESSURE NORMAL TO THE SURFACE

  // auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg->vertexNormals);
  // externalPressure_e = externalPressureMagnitude *
  // vertexAngleNormal_e.row(P.ptInd);

  // c. ALTERNATIVELY, PRESSURE BASED ON INITIAL GEOMETRY + ALONG A FIXED
  // DIRECTION, E.G. NEGATIVE Z DIRECTION

  // initialize/update the external pressure magnitude distribution
  gaussianDistribution(externalPressureMagnitude,
                       geodesicDistanceFromPtInd.raw(),
                       geodesicDistanceFromPtInd.raw().maxCoeff() / P.conc);
  externalPressureMagnitude *= P.Kf;

  Eigen::Matrix<double, 1, 3> zDir;
  zDir << 0.0, 0.0, -1.0;
  // externalPressure_e = -externalPressureMagnitude * zDir *
  //                      (vpg->inputVertexPositions[theVertex].z - P.height);
  F.externalForce.raw() = externalPressureMagnitude;
  F.externalForce.raw() *= vpg->vertexLumpedMassMatrix;

  return F.externalForce.raw();
}

EigenVectorX1d System::computeChemicalPotential() {
  gcs::VertexData<double> dH0dphi(*mesh, 0);
  gcs::VertexData<double> dKbdphi(*mesh, 0);
  auto meanCurvDiff = (vpg->vertexMeanCurvatures.raw().array() /
                       vpg->vertexDualAreas.raw().array()) -
                      H0.raw().array();

  if (P.relation == "linear") {
    dH0dphi.fill(P.H0c);
    dKbdphi.fill(P.Kbc);
  } else if (P.relation == "hill") {
    EigenVectorX1d proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    dH0dphi.raw() =
        (2 * P.H0c * proteinDensity.raw().array() /
         ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
            .matrix();
    dKbdphi.raw() =
        (2 * P.Kbc * proteinDensity.raw().array() /
         ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
            .matrix();
  }

  F.adsorptionPotential.raw() =
      F.maskProtein(-P.epsilon * vpg->vertexDualAreas.raw().array());
  F.bendingPotential.raw() = F.maskProtein(
      -vpg->vertexDualAreas.raw().array() *
      (meanCurvDiff * meanCurvDiff * dKbdphi.raw().array() -
       2 * Kb.raw().array() * meanCurvDiff * dH0dphi.raw().array()));
  F.diffusionPotential.raw() =
      F.maskProtein(-P.eta * vpg->cotanLaplacian * proteinDensity.raw());
  F.interiorPenaltyPotential.raw() =
      F.maskProtein(P.lambdaPhi * (1 / proteinDensity.raw().array() -
                                   1 / (1 - proteinDensity.raw().array())));
  F.chemicalPotential.raw() =
      F.adsorptionPotential.raw() + F.bendingPotential.raw() +
      F.diffusionPotential.raw() + F.interiorPenaltyPotential.raw();

  // F.chemicalPotential.raw().array() =
  //     -vpg->vertexDualAreas.raw().array() *
  //     (P.epsilon - 2 * Kb.raw().array() * meanCurvDiff *
  //     dH0dphi.raw().array() +
  //      meanCurvDiff * meanCurvDiff * dKbdphi.raw().array());
  // F.chemicalPotential.raw().array() +=
  //     P.lambdaPhi * (1 / proteinDensity.raw().array() -
  //                    1 / (1 - proteinDensity.raw().array()));

  return F.chemicalPotential.raw();
}

std::tuple<EigenVectorX3dr, EigenVectorX3dr>
System::computeDPDForces(double dt) {

  auto dampingForce_e = EigenMap<double, 3>(F.dampingForce);
  auto stochasticForce_e = EigenMap<double, 3>(F.stochasticForce);

  // Reset forces to zero
  dampingForce_e.setZero();
  stochasticForce_e.setZero();

  // alias positions
  const auto &pos = vpg->inputVertexPositions;

  // std::default_random_engine random_generator;
  // gcs::EdgeData<double> random_var(mesh);
  double sigma =
      sqrt(2 * P.gamma * mem3dg::constants::kBoltzmann * P.temp / dt);
  std::normal_distribution<double> normal_dist(0, sigma);

  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    gcs::Vertex v1 = he.vertex();
    gcs::Vertex v2 = he.next().vertex();

    gc::Vector3 dVel12 = vel[v1] - vel[v2];
    gc::Vector3 dPos12_n = (pos[v1] - pos[v2]).normalize();

    if (P.gamma != 0) {
      gc::Vector3 df = P.gamma * (gc::dot(dVel12, dPos12_n) * dPos12_n);
      F.dampingForce[v1] -= df;
      F.dampingForce[v2] += df;
    }

    if (sigma != 0) {
      double noise = normal_dist(rng);
      F.stochasticForce[v1] += noise * dPos12_n;
      F.stochasticForce[v2] -= noise * dPos12_n;
    }

    // gc::Vector3 dVel21 = vel[v2] - vel[v1];
    // gc::Vector3 dPos21_n = (pos[v2] - pos[v1]).normalize();

    // std::cout << -gamma * (gc::dot(dVel12, dPos12_n) * dPos12_n)
    //           << " == " << -gamma * (gc::dot(-dVel12, -dPos12_n) * -dPos12_n)
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
  double q1 = quantities[he.next().next().vertex()];
  double q2 = quantities[he.vertex()];
  double q3 = quantities[he.next().vertex()];

  if (q1 == q2 && q1 == q3) {
    return gc::Vector3({0, 0, 0});
  } else {
    // Edge and normal vector
    gc::Vector3 n = vpg->faceNormals[he.face()];
    gc::Vector3 e1 = vecFromHalfedge(he, *vpg);
    gc::Vector3 e2 = vecFromHalfedge(he.next(), *vpg);
    gc::Vector3 e3 = vecFromHalfedge(he.next().next(), *vpg);

    // exterior angle of triangles (angles formed by e_perp)
    double theta3 = gc::angle(e1, e2);
    double theta1 = gc::angle(e2, e3);
    double theta2 = gc::angle(e3, e1);

    // gradient of edge length wrt he.vertex()
    gc::Vector3 grad_e1norm = -e1.normalize();
    gc::Vector3 grad_e3norm = e3.normalize();

    // gradient of exterior angle wrt he.vertex()
    gc::Vector3 grad_theta3 = gc::cross(n, e1).normalize() / gc::norm(e1);
    gc::Vector3 grad_theta1 = gc::cross(n, e3).normalize() / gc::norm(e3);
    gc::Vector3 grad_theta2 = -(grad_theta3 + grad_theta1);

    // chain rule
    gc::Vector3 grad_costheta3 = -sin(theta3) * grad_theta3;
    gc::Vector3 grad_costheta2 = -sin(theta2) * grad_theta2;
    gc::Vector3 grad_costheta1 = -sin(theta1) * grad_theta1;

    // g = q1 * e1_perp +  q2 * e2_perp +  q2 * e2_perp
    // gradient of |g|^2
    return 2 * q1 * q1 * gc::norm(e1) * grad_e1norm +
           2 * q3 * q3 * gc::norm(e3) * grad_e3norm +
           2 * q1 * q2 * gc::norm(e2) *
               (grad_e1norm * cos(theta3) + gc::norm(e1) * grad_costheta3) +
           2 * q2 * q3 * gc::norm(e2) *
               (grad_e3norm * cos(theta1) + gc::norm(e3) * grad_costheta1) +
           2 * q1 * q3 *
               (grad_e1norm * gc::norm(e3) * cos(theta2) +
                gc::norm(e1) * grad_e3norm * cos(theta2) +
                gc::norm(e1) * gc::norm(e3) * grad_costheta2);
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
  F.bendingForceVec.fill({0, 0, 0});
  F.capillaryForceVec.fill({0, 0, 0});
  F.osmoticForceVec.fill({0, 0, 0});
  F.lineCapillaryForceVec.fill({0, 0, 0});
  F.adsorptionForceVec.fill({0, 0, 0});

  F.bendingForce.raw().setZero();
  F.capillaryForce.raw().setZero();
  F.lineCapillaryForce.raw().setZero();
  F.externalForce.raw().setZero();
  F.osmoticForce.raw().setZero();

  F.chemicalPotential.raw().setZero();
  F.diffusionPotential.raw().setZero();
  F.bendingPotential.raw().setZero();
  F.adsorptionPotential.raw().setZero();
  F.interiorPenaltyPotential.raw().setZero();

  if (O.isShapeVariation) {
    computeVectorForces();
    if (P.Kf != 0) {
      computeExternalForce();
    }
    F.mechanicalForceVec = F.osmoticForceVec + F.capillaryForceVec +
                           F.bendingForceVec + F.lineCapillaryForceVec +
                           F.adsorptionForceVec + F.addNormal(F.externalForce);
    F.mechanicalForce = F.ontoNormal(F.mechanicalForceVec);
  }

  if (O.isProteinVariation) {
    computeChemicalPotential();
  }

  // computeBendingForce();
  // if (P.Kv != 0) {
  //   computeOsmoticForce();
  // }
  // if (P.Ksg != 0) {
  //   computeCapillaryForce();
  // }
  // if (P.eta != 0) {
  //   computeLineCapillaryForce();
  // }

  // compute the mechanical error norm
  mechErrorNorm =
      O.isShapeVariation ? computeNorm(toMatrix(F.mechanicalForceVec)) : 0;

  // compute the chemical error norm
  chemErrorNorm =
      O.isProteinVariation ? computeNorm(F.chemicalPotential.raw()) : 0;
}

} // namespace solver
} // namespace mem3dg
