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

#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/solver/system.h"
#include "mem3dg/solver/meshops.h"
#include <Eigen/Core>

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::getBendingPressure() {

  // map the MeshData to eigen matrix XXX_e
  auto bendingPressure_e = gc::EigenMap<double, 3>(bendingPressure);
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
  auto positions = gc::EigenMap<double, 3>(vpg.inputVertexPositions);

  // Alias
  std::size_t n_vertices = (mesh.nVertices());

  /// A. non-optimized version
  // calculate the Laplacian of mean curvature H
  Eigen::Matrix<double, Eigen::Dynamic, 1> lap_H = M_inv * L * (H - H0);

  // initialize and calculate intermediary result scalerTerms
  Eigen::Matrix<double, Eigen::Dynamic, 1> scalerTerms =
      rowwiseProduct(H, H) + rowwiseProduct(H, H0) -
      M_inv * vpg.vertexGaussianCurvatures.raw();
  Eigen::Matrix<double, Eigen::Dynamic, 1> zeroMatrix;
  zeroMatrix.resize(n_vertices, 1);
  zeroMatrix.setZero();
  scalerTerms = scalerTerms.array().max(zeroMatrix.array());

  // initialize and calculate intermediary result productTerms
  Eigen::Matrix<double, Eigen::Dynamic, 1> productTerms;
  productTerms.resize(n_vertices, 1);
  productTerms = 2.0 * rowwiseProduct(scalerTerms, H - H0);

  // calculate bendingForce
  bendingPressure_e =
      -2.0 * P.Kb * rowwiseScaling(productTerms + lap_H, vertexAngleNormal_e);

  // /// B. optimized version
  // // calculate the Laplacian of mean curvature H
  // Eigen::Matrix<double, Eigen::Dynamic, 1> lap_H_integrated = L * (H - H0);

  // // initialize and calculate intermediary result scalarTerms_integrated
  // Eigen::Matrix<double, Eigen::Dynamic, 1> H_integrated = M * H;
  // Eigen::Matrix<double, Eigen::Dynamic, 1> scalarTerms_integrated =
  //     M * rowwiseProduct(M_inv * H_integrated, M_inv * H_integrated) +
  //     rowwiseProduct(H_integrated, H0) - vpg.vertexGaussianCurvatures.raw();
  // Eigen::Matrix<double, Eigen::Dynamic, 1> zeroMatrix;
  // zeroMatrix.resize(n_vertices, 1);
  // zeroMatrix.setZero();
  // scalarTerms_integrated =
  //     scalarTerms_integrated.array().max(zeroMatrix.array());

  // // initialize and calculate intermediary result productTerms_integrated
  // Eigen::Matrix<double, Eigen::Dynamic, 1> productTerms_integrated;
  // productTerms_integrated.resize(n_vertices, 1);
  // productTerms_integrated =
  //     2.0 * rowwiseProduct(scalarTerms_integrated, H - H0);

  // bendingPressure_e =
  //     -2.0 * P.Kb *
  //     rowwiseScaling(M_inv * (productTerms_integrated + lap_H_integrated),
  //                    vertexAngleNormal_e);
}

void System::getCapillaryPressure() {

  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
  auto capillaryPressure_e = gc::EigenMap<double, 3>(capillaryPressure);

  /// Geometric implementation
  if (mesh.hasBoundary()) {
    /// CAPILLARY PRESSURE OF PATCH
    capillaryPressure_e = rowwiseScaling(-P.Ksg * 2.0 * H, vertexAngleNormal_e);
  } else {
    /// CAPILLARY PRESSURE OF VESICLE
    capillaryPressure_e = rowwiseScaling(
        -(P.Ksg * (surfaceArea - targetSurfaceArea) / targetSurfaceArea +
          P.lambdaSG) *
            2.0 * H,
        vertexAngleNormal_e);
  }

  // /// Nongeometric implementation
  // for (gcs::Vertex v : mesh.vertices()) {
  //   gc::Vector3 globalForce{0.0, 0.0, 0.0};
  //   for (gcs::Halfedge he : v.outgoingHalfedges()) {
  //     gc::Vector3 base_vec = vecFromHalfedge(he.next(), vpg);
  //     gc::Vector3 localAreaGradient =
  //         -gc::cross(base_vec, vpg.faceNormals[he.face()]);
  //     assert((gc::dot(localAreaGradient, vecFromHalfedge(he, vpg))) < 0);
  //     if (P.Ksg != 0) {
  //       capillaryPressure[v] += -P.Ksg * localAreaGradient *
  //                               (surfaceArea - targetSurfaceArea) /
  //                               targetSurfaceArea;
  //     }
  //   }
  //   capillaryPressure[v] /= vpg.vertexDualAreas[v];
  // }
}

void System::getInsidePressure() {
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
  auto insidePressure_e = gc::EigenMap<double, 3>(insidePressure);

  /// Geometric implementation
  if (mesh.hasBoundary()) {
    /// Inside excess pressure of patch
    insidePressure_e = P.Kv * vertexAngleNormal_e;
  } else {
    /// Inside excess pressure of vesicle
    insidePressure_e =
        -(P.Kv * (volume - refVolume * P.Vt) / (refVolume * P.Vt) + P.lambdaV) *
        vertexAngleNormal_e;
  }

  // /// Nongeometric implementation
  // for (gcs::Vertex v : mesh.vertices()) {
  //   for (gcs::Halfedge he : v.outgoingHalfedges()) {
  //     gc::Vector3 p1 = vpg.inputVertexPositions[he.next().vertex()];
  //     gc::Vector3 p2 = vpg.inputVertexPositions[he.next().next().vertex()];
  //     gc::Vector3 dVdx = 0.5 * gc::cross(p1, p2) / 3.0;
  //     insidePressure[v] +=
  //         -P.Kv * (volume - refVolume * P.Vt) / (refVolume * P.Vt) * dVdx;
  //   }
  // }
}

void System::getLineTensionPressure() {
  gcs::HalfedgeData<gc::Vector2> halfedgeVectorsInVertex(mesh);
  gcs::VertexData<std::array<gc::Vector3, 2>> vertexTangentBasis(mesh);

  for (gcs::Vertex v : mesh.vertices()) {
    // Calculate interfacial tension
    if ((H0[v.getIndex()] > (0.1 * P.H0)) &&
        (H0[v.getIndex()] < (0.9 * P.H0)) && (H[v.getIndex()] != 0)) {

      // Calculate halfedgeVectorsInVertex
      double coordSum = 0.0;
      gcs::Halfedge firstHe = v.halfedge();
      gcs::Halfedge currHe = firstHe;
      do {
        halfedgeVectorsInVertex[currHe] =
            gc::Vector2::fromAngle(coordSum) * vpg.edgeLengths[currHe.edge()];
        coordSum += vpg.cornerScaledAngles[currHe.corner()];
        if (!currHe.isInterior())
          break;
        currHe = currHe.next().next().twin();
      } while (currHe != firstHe);

      // Initialize principal direction, gradient and tagentBasis
      gc::Vector2 principalDirection1{0.0, 0.0};
      gc::Vector3 gradient{0.0, 0.0, 0.0};
      gc::Vector3 basisXSum = gc::Vector3::zero();

      // Calculate principal direction, gradient and tangentBasis
      for (gcs::Halfedge he : v.outgoingHalfedges()) {

        // Calculate dihedral angle alpha
        gcs::Edge e = he.edge();
        if (e.isBoundary())
          continue;
        if (!e.isManifold()) {
          continue;
        }
        gc::Vector3 N1 = vpg.faceNormals[e.halfedge().face()];
        gc::Vector3 N2 = vpg.faceNormals[e.halfedge().sibling().face()];
        gc::Vector3 pTail = vpg.vertexPositions[e.halfedge().vertex()];
        gc::Vector3 pTip = vpg.vertexPositions[e.halfedge().next().vertex()];
        gc::Vector3 edgeDir = gc::unit(pTip - pTail);
        double alpha = std::atan2(dot(edgeDir, cross(N1, N2)), dot(N1, N2));

        // calculate principal direction
        double len = vpg.edgeLengths[he.edge()];
        gc::Vector2 vec = halfedgeVectorsInVertex[he];
        principalDirection1 += -vec * vec / len * std::abs(alpha);

        // calculate gradient
        gradient +=
            vecFromHalfedge(he, vpg).normalize() *
            (H0[he.next().vertex().getIndex()] - H0[he.vertex().getIndex()]) /
            vpg.edgeLengths[he.edge()];

        // calculate tangent basis
        gc::Vector3 eVec = vpg.vertexPositions[he.next().vertex()] -
                           vpg.vertexPositions[he.vertex()];
        eVec = eVec.removeComponent(vpg.vertexNormals[v]);
        double angle = halfedgeVectorsInVertex[he].arg();
        gc::Vector3 eVecX = eVec.rotateAround(vpg.vertexNormals[v], -angle);
        basisXSum += eVecX;

      }

      // post-process gradient and vertex principal direction
      gradient.normalize();
      principalDirection1 /= 4.0;
      gc::Vector3 basisX = unit(basisXSum);
      gc::Vector3 basisY = cross(vpg.vertexNormals[v], basisX);
      vertexTangentBasis[v] = {{basisX, basisY}};

      // Find angle between tangent & principal direction
      gc::Vector3 tangentVector =
          gc::cross(gradient, vpg.vertexNormals[v]).normalize();
      gc::Vector3 PD1InWorldCoords =
          vertexTangentBasis[v][0] * principalDirection1.x +
          vertexTangentBasis[v][1] * principalDirection1.y;
      double cosT = gc::dot(tangentVector, PD1InWorldCoords.normalize());

      // Deduce normal curvature
      double K1 =
          (2 * H[v.getIndex()] + sqrt(principalDirection1.norm())) * 0.5;
      double K2 =
          (2 * H[v.getIndex()] - sqrt(principalDirection1.norm())) * 0.5;
      lineTensionPressure[v] = -P.eta * vpg.vertexNormals[v] *
                               (cosT * cosT * (K1 - K2) + K2) * P.sharpness;
    }
  }

  // /// If requireVertexPrincipalCurvatureDirections and requireVertexTangentBasis
  // for (gcs::Vertex v : mesh.vertices()) {
  //   // Calculate interfacial tension
  //   if ((H0[v.getIndex()] > (0.1 * P.H0)) &&
  //       (H0[v.getIndex()] < (0.9 * P.H0)) && (H[v.getIndex()] != 0)) {
  //     gc::Vector3 gradient{0.0, 0.0, 0.0};
  //     // Calculate gradient of spon curv
  //     for (gcs::Halfedge he : v.outgoingHalfedges()) {
  //       gradient +=
  //           vecFromHalfedge(he, vpg).normalize() *
  //           (H0[he.next().vertex().getIndex()] - H0[he.vertex().getIndex()]) /
  //           vpg.edgeLengths[he.edge()];
  //     }
  //     gradient.normalize();
  //     // Find angle between tangent & principal direction
  //     gc::Vector3 tangentVector =
  //         gc::cross(gradient, vpg.vertexNormals[v]).normalize();
  //     gc::Vector2 principalDirection1 =
  //         vpg.vertexPrincipalCurvatureDirections[v];
  //     gc::Vector3 PD1InWorldCoords =
  //         vpg.vertexTangentBasis[v][0] * principalDirection1.x +
  //         vpg.vertexTangentBasis[v][1] * principalDirection1.y;
  //     double cosT = gc::dot(tangentVector, PD1InWorldCoords.normalize());
  //     // Deduce normal curvature
  //     double K1 =
  //         (2 * H[v.getIndex()] + sqrt(principalDirection1.norm())) * 0.5;
  //     double K2 =
  //         (2 * H[v.getIndex()] - sqrt(principalDirection1.norm())) * 0.5;
  //     lineTensionPressure[v] = -P.eta * vpg.vertexNormals[v] *
  //                              (cosT * cosT * (K1 - K2) + K2) * P.sharpness;
  //   }
  // }

}

void System::getExternalPressure() {

  auto externalPressure_e = gc::EigenMap<double, 3>(externalPressure);
  Eigen::Matrix<double, Eigen::Dynamic, 1> externalPressureMagnitude;

  if (P.Kf != 0) {

    // a. FIND OUT THE CURRENT EXTERNAL PRESSURE MAGNITUDE BASED ON CURRENT
    // GEOMETRY

    // auto &dist_e = heatMethodDistance(vpg, mesh.vertex(P.ptInd)).raw();
    // double stdDev = dist_e.maxCoeff() / P.conc;
    // externalPressureMagnitude =
    //    P.Kf / (stdDev * pow(M_PI * 2, 0.5)) *
    //    (-dist_e.array() * dist_e.array() / (2 * stdDev * stdDev)).exp();

    // b. APPLY EXTERNAL PRESSURE NORMAL TO THE SURFACE

    // auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
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
    externalPressure_e =
        -externalPressureMagnitude * zDir *
        (vpg.inputVertexPositions[mesh.vertex(ptInd)].z - P.height);
  }
}

void System::getChemicalPotential() {

  Eigen::Matrix<double, Eigen::Dynamic, 1> proteinDensitySq =
      (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();

  Eigen::Matrix<double, Eigen::Dynamic, 1> dH0dphi =
      (2 * P.H0 * proteinDensity.raw().array() /
       ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
          .matrix();

  chemicalPotential.raw() =
      (P.epsilon - (2 * P.Kb * (H - H0)).array() * dH0dphi.array()).matrix();
}

} // namespace ddgsolver