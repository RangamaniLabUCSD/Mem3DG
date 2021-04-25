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
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"
#include "geometrycentral/utilities/vector3.h"
#include "mem3dg/solver/constants.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"
#include <Eigen/Core>
#include <cmath>

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

gcs::EdgeData<double>
System::computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg) const {
  gcs::EdgeData<double> LCR(*mesh);
  for (gcs::Edge e : mesh->edges()) {
    gcs::Edge lj = e.halfedge().next().edge();
    gcs::Edge ki = e.halfedge().twin().next().edge();
    gcs::Edge il = e.halfedge().next().next().edge();
    gcs::Edge jk = e.halfedge().twin().next().next().edge();
    LCR[e] = vpg.edgeLengths[il] * vpg.edgeLengths[jk] / vpg.edgeLengths[ki] /
             vpg.edgeLengths[lj];
  }
  return LCR;
}

double System::computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                                       gcs::Edge &e) const {

  gcs::Edge lj = e.halfedge().next().edge();
  gcs::Edge ki = e.halfedge().twin().next().edge();
  gcs::Edge il = e.halfedge().next().next().edge();
  gcs::Edge jk = e.halfedge().twin().next().next().edge();
  return vpg.edgeLengths[il] * vpg.edgeLengths[jk] / vpg.edgeLengths[ki] /
         vpg.edgeLengths[lj];
}

double System::computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                                       gcs::Edge &&e) const {

  gcs::Edge lj = e.halfedge().next().edge();
  gcs::Edge ki = e.halfedge().twin().next().edge();
  gcs::Edge il = e.halfedge().next().next().edge();
  gcs::Edge jk = e.halfedge().twin().next().next().edge();
  return vpg.edgeLengths[il] * vpg.edgeLengths[jk] / vpg.edgeLengths[ki] /
         vpg.edgeLengths[lj];
}

void System::computeRegularizationForce() {
  for (gcs::Vertex v : mesh->vertices()) {
    if (!v.isBoundary()) {
      for (gcs::Halfedge he : v.outgoingHalfedges()) {

        // Conformal regularization
        if (P.Kst != 0 && !he.edge().isBoundary()) {
          gcs::Halfedge jl = he.next();
          gcs::Halfedge li = jl.next();
          gcs::Halfedge ik = he.twin().next();
          gcs::Halfedge kj = ik.next();

          gc::Vector3 grad_li = vecFromHalfedge(li, *vpg).normalize();
          gc::Vector3 grad_ik = vecFromHalfedge(ik.twin(), *vpg).normalize();
          F.regularizationForce[v] +=
              -P.Kst *
              (computeLengthCrossRatio(*vpg, he.edge()) -
               targetLcrs[he.edge()]) /
              targetLcrs[he.edge()] *
              (vpg->edgeLengths[kj.edge()] / vpg->edgeLengths[jl.edge()]) *
              (grad_li * vpg->edgeLengths[ik.edge()] -
               grad_ik * vpg->edgeLengths[li.edge()]) /
              vpg->edgeLengths[ik.edge()] / vpg->edgeLengths[ik.edge()];
        }

        // Local area regularization
        if (P.Ksl != 0 && he.isInterior()) {
          gcs::Halfedge base_he = he.next();
          gc::Vector3 base_vec = vecFromHalfedge(base_he, *vpg);
          gc::Vector3 localAreaGradient =
              -gc::cross(base_vec, vpg->faceNormals[he.face()]);
          auto &referenceArea = (v.isBoundary() ? refFaceAreas[base_he.face()]
                                                : meanTargetFaceArea);
          F.regularizationForce[v] +=
              -P.Ksl * localAreaGradient *
              (vpg->faceAreas[base_he.face()] - referenceArea);
        }

        // local edge regularization
        if (P.Kse != 0) {
          gc::Vector3 edgeGradient = -vecFromHalfedge(he, *vpg).normalize();
          auto &referenceLength = (v.isBoundary() ? refEdgeLengths[he.edge()]
                                                  : meanTargetEdgeLength);
          F.regularizationForce[v] +=
              -P.Kse * edgeGradient *
              (vpg->edgeLengths[he.edge()] - referenceLength);
        }
      }
    }
  }

  // post processing regularization force
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg->vertexNormals);
  auto regularizationForce_e = gc::EigenMap<double, 3>(F.regularizationForce);

  // remove the normal component
  regularizationForce_e -= rowwiseScaling(
      rowwiseDotProduct(regularizationForce_e, vertexAngleNormal_e),
      vertexAngleNormal_e);

  // // moving boundary
  // for (gcs::Vertex v : mesh->vertices()) {
  //   if (v.isBoundary()) {
  //     F.regularizationForce[v].z = 0;
  //   }
  // }

  // / Patch regularization
  // the cubic penalty is for regularizing the mesh,
  // need better physical interpretation or alternative method
  // if (P.Kse != 0) {
  //   double strain =
  //       (vpg->edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
  //       targetEdgeLengths[he.edge()];
  //   F.regularizationForce[v] +=
  //      -P.Kse * edgeGradient * strain * strain * strain;
  // }

  // // remove the masked components
  // regularizationForce_e =
  //     rowwiseScaling(mask.raw().cast<double>(), regularizationForce_e);
  // // moving boundary
  // for (gcs::Vertex v : mesh->vertices()) {
  //   if (!mask[v]) {
  //     F.regularizationForce[v].z = 0;
  //     // boundary tension, mostly likely not necessary
  //     if (v.isBoundary()) {
  //       double boundaryEdgeLength = 0;
  //       for (gcs::Edge e : v.adjacentEdges()) {
  //         if (e.isBoundary()) {
  //           boundaryEdgeLength += vpg->edgeLength(e);
  //         }
  //       }
  //       boundaryEdgeLength /= 2;
  //       double scaling;
  //       if (F.regularizationForce[v].norm() > 1e-15) {
  //         scaling = 1 - abs(P.Ksg * boundaryEdgeLength /
  //                           (F.regularizationForce[v].norm()));
  //       }
  //       F.regularizationForce[v] *= scaling;
  //     }
  //   }
  // }
}

void System::vertexShift() {
  for (gcs::Vertex v : mesh->vertices()) {
    if (mask[v]) {
      gc::Vector3 baryCenter{0.0, 0.0, 0.0};
      double n_vAdj = 0.0;
      for (gcs::Vertex vAdj : v.adjacentVertices()) {
        baryCenter += vpg->inputVertexPositions[vAdj];
        n_vAdj += 1.0;
      }
      baryCenter /= n_vAdj;
      for (gcs::Halfedge he : v.outgoingHalfedges()) {
        gcs::Halfedge base_he = he.next();
        vpg->inputVertexPositions[v] =
            baryCenter - gc::dot(vpg->vertexNormals[v],
                                 baryCenter - vpg->inputVertexPositions[v]) *
                             vpg->vertexNormals[v];
      }
    }
  }
}

bool System::edgeFlip() {
  bool isFlipped = false;

  // flip edge if not delauney
  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    // if (mask[he.vertex()] || mask[he.twin().vertex()]) {
    if (!he.edge().isBoundary()) {
      bool nonDelaunay =
          (vpg->cornerAngle(he.next().next().corner()) +
           vpg->cornerAngle(he.twin().next().next().corner())) > constants::PI;
      bool flat = abs(vpg->edgeDihedralAngle(he.edge())) < (constants::PI / 36);
      if (nonDelaunay && flat) {
        auto success = mesh->flip(he.edge());
        isFlipped = true;

        auto maskAllNeighboring = [](gcs::VertexData<bool> &smoothingMask,
                                     const gcs::Vertex v) {
          smoothingMask[v] = true;
          for (gc::Vertex nv : v.adjacentVertices()) {
            smoothingMask[nv] = true;
          }
        };
        maskAllNeighboring(smoothingMask, he.tailVertex());
        maskAllNeighboring(smoothingMask, he.tipVertex());
      }
    }
  }

  if (isFlipped)
    mesh->compress();

  return isFlipped;
}

bool System::growMesh() {
  bool isGrown = false;
  // expand the mesh when area is too large
  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    // if (mask[he.vertex()] || mask[he.twin().vertex()]) {
    if (!he.edge().isBoundary()) {
      // curvature based remeshing:
      // https://www.irit.fr/recherches/VORTEX/publications/rendu-geometrie/EGshort2013_Dunyach_et_al.pdf
      double k1 = (abs(vpg->vertexMaxPrincipalCurvature(he.tipVertex())) >
                   abs(vpg->vertexMinPrincipalCurvature(he.tipVertex())))
                      ? abs(vpg->vertexMaxPrincipalCurvature(he.tipVertex()))
                      : abs(vpg->vertexMinPrincipalCurvature(he.tipVertex()));
      double k2 = (abs(vpg->vertexMaxPrincipalCurvature(he.tailVertex())) >
                   abs(vpg->vertexMinPrincipalCurvature(he.tailVertex())))
                      ? abs(vpg->vertexMaxPrincipalCurvature(he.tailVertex()))
                      : abs(vpg->vertexMinPrincipalCurvature(he.tailVertex()));
      double L = std::sqrt(6 * P.curvTol / ((k1 > k2) ? k1 : k2) -
                           3 * P.curvTol * P.curvTol);
      bool is2Large =
          (vpg->faceArea(he.face()) + vpg->faceArea(he.twin().face())) >
          (4 * meanTargetFaceArea);
      bool is2Curved = vpg->edgeLength(he.edge()) > (2 * L);
      bool nonDelaunay =
          (vpg->cornerAngle(he.next().next().corner()) +
           vpg->cornerAngle(he.twin().next().next().corner())) > constants::PI;
      bool flat = abs(vpg->edgeDihedralAngle(he.edge())) < (constants::PI / 36);

      if (is2Large || is2Curved || (nonDelaunay && !flat)) {
        // || abs(H0[he.tailVertex()] - H0[he.tipVertex()]) > 0.2) {

        // split the edge
        const auto &vertex1 = he.tipVertex(), &vertex2 = he.tailVertex();
        const auto &newVertex = mesh->splitEdgeTriangular(he.edge()).vertex();

        // update quantities
        averageData(vpg->inputVertexPositions, vertex1, vertex2, newVertex);
        // Note: think about conservation of energy, momentum and angular
        // momentum
        averageData(vel, vertex1, vertex2, newVertex);
        averageData(geodesicDistanceFromPtInd, vertex1, vertex2, newVertex);
        thePointTracker[newVertex] = false;
        if (O.isProtein)
          averageData(proteinDensity, vertex1, vertex2, newVertex);

        // smoothing mask
        auto maskAllNeighboring = [](gcs::VertexData<bool> &smoothingMask,
                                     const gcs::Vertex v) {
          smoothingMask[v] = true;
          for (gc::Vertex nv : v.adjacentVertices()) {
            smoothingMask[nv] = true;
          }
        };
        maskAllNeighboring(smoothingMask, newVertex);
        // smoothingMask[newVertex] = true;

        isGrown = true;
      }
    }
  }

  // shrink the mesh is vertex is too close
  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    // if (mask[he.vertex()] || mask[he.twin().vertex()]) {
    if (!he.edge().isBoundary()) {
      bool is2Skinny = (vpg->cornerAngle(he.next().next().corner()) +
                        vpg->cornerAngle(he.twin().next().next().corner())) <
                       constants::PI / 3;
      if (is2Skinny) {

        // alias the neighboring vertices
        const auto &vertex1 = he.tipVertex(), &vertex2 = he.tailVertex();

        // precached pre-mutation values or flag
        gc::Vector3 collapsedPosition =
            vertex1.isBoundary()   ? vpg->inputVertexPositions[vertex1]
            : vertex2.isBoundary() ? vpg->inputVertexPositions[vertex2]
                                   : (vpg->inputVertexPositions[vertex1] +
                                      vpg->inputVertexPositions[vertex2]) /
                                         2;
        bool isThePoint = thePointTracker[vertex1] || thePointTracker[vertex2];

        // collapse the edge
        auto newVertex = mesh->collapseEdgeTriangular(he.edge());

        // update quantities
        vpg->inputVertexPositions[newVertex] = collapsedPosition;
        thePointTracker[newVertex] = isThePoint;
        // Note: think about conservation of energy, momentum and angular
        // momentum
        averageData(vel, vertex1, vertex2, newVertex);
        averageData(geodesicDistanceFromPtInd, vertex1, vertex2, newVertex);
        if (O.isProtein)
          averageData(proteinDensity, vertex1, vertex2, newVertex);

        // smoothing mask
        auto maskAllNeighboring = [](gcs::VertexData<bool> &smoothingMask,
                                     const gcs::Vertex v) {
          smoothingMask[v] = true;
          for (gc::Vertex nv : v.adjacentVertices()) {
            smoothingMask[nv] = true;
          }
        };
        maskAllNeighboring(smoothingMask, newVertex);

        isGrown = true;
        // std::cout << "in shirnk sum: "
        //           << thePointTracker.raw().cast<double>().sum() <<
        //           std::endl;
      }
    }
  }

  if (isGrown)
    mesh->compress();
  return isGrown;
}

void System::processMesh() {

  bool isGrown = false, isFlipped = false;
  smoothingMask.fill(false);

  // vertex shift for regularization
  if (O.isVertexShift) {
    vertexShift();
  }

  // split edge and collapse edge
  if (O.isGrowMesh) {
    isGrown = growMesh();
  }

  // linear edge flip for non-Delauney triangles
  if (O.isEdgeFlip) {
    isFlipped = edgeFlip();
  }

  // regularization
  if ((P.Kse != 0) || (P.Ksl != 0) || (P.Kst != 0)) {
    computeRegularizationForce();
    vpg->inputVertexPositions.raw() += F.regularizationForce.raw();
  }

  // globally update quantities
  if (isGrown || isFlipped) {
    // globalSmoothing(smoothingMask);
    globalUpdateAfterMutation();
  }
}

void System::globalSmoothing(gcs::VertexData<bool> &smoothingMask, double tol,
                             double stepSize) {
  EigenVectorX1D gradient;
  double pastGradNorm = 1e10;
  double gradNorm;
  do {
    vpg->refreshQuantities();
    computeBendingForce();
    auto pos_e = gc::EigenMap<double, 3>(vpg->inputVertexPositions);
    auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg->vertexNormals);
    gradient = (smoothingMask.raw().cast<double>()).array() *
               F.bendingForce.raw().array();
    gradNorm =
        gradient.cwiseAbs().sum() / smoothingMask.raw().cast<int>().sum();
    if (gradNorm > pastGradNorm) {
      stepSize /= 2;
      // std::cout << "WARNING: globalSmoothing: stepSize too large, cut in
      // half!"
      //           << std::endl;
    }
    pos_e.array() +=
        rowwiseScaling(gradient, vertexAngleNormal_e).array() * stepSize;
    pastGradNorm = gradNorm;
    // std::cout << "gradient:  " << gradNorm << std::endl;
  } while (gradNorm > tol);
}

void System::localSmoothing(const gcs::Vertex &v, std::size_t num,
                            double stepSize) {
  std::size_t count = 0;
  while (count < num) {
    gc::Vector3 vertexNormal{0, 0, 0};
    for (gcs::Corner c : v.adjacentCorners()) {
      vertexNormal += vpg->cornerAngle(c) * vpg->faceNormal(c.face());
    }
    vertexNormal.normalize();
    double localLapH = 0;
    double H_center = vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v);
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      localLapH += vpg->edgeCotanWeight(he.edge()) *
                   (H_center - vpg->vertexMeanCurvature(he.tipVertex()) /
                                   vpg->vertexDualArea(he.tipVertex()));
    }
    vpg->inputVertexPositions[v] -=
        stepSize * vpg->vertexDualArea(v) * localLapH * vertexNormal;
    count++;
  }
}

void System::localSmoothing(const gcs::Halfedge &he, std::size_t num,
                            double stepSize) {
  std::size_t count = 0;
  while (count < num) {
    gc::Vector3 vertexNormal1{0, 0, 0};
    gc::Vector3 vertexNormal2{0, 0, 0};
    double localLapH1 = 0;
    double localLapH2 = 0;

    auto v = he.tailVertex();
    for (gcs::Corner c : v.adjacentCorners()) {
      vertexNormal1 += vpg->cornerAngle(c) * vpg->faceNormal(c.face());
    }
    vertexNormal1.normalize();
    double H_center = vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v);
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      localLapH1 += vpg->edgeCotanWeight(he.edge()) *
                    (H_center - vpg->vertexMeanCurvature(he.tipVertex()) /
                                    vpg->vertexDualArea(he.tipVertex()));
    }

    v = he.tipVertex();
    for (gcs::Corner c : v.adjacentCorners()) {
      vertexNormal2 += vpg->cornerAngle(c) * vpg->faceNormal(c.face());
    }
    vertexNormal2.normalize();
    H_center = vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v);
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      localLapH2 += vpg->edgeCotanWeight(he.edge()) *
                    (H_center - vpg->vertexMeanCurvature(he.tipVertex()) /
                                    vpg->vertexDualArea(he.tipVertex()));
    }

    vpg->inputVertexPositions[he.tailVertex()] -=
        stepSize * vpg->vertexDualArea(he.tailVertex()) * localLapH1 *
        vertexNormal1;
    vpg->inputVertexPositions[he.tipVertex()] -=
        stepSize * vpg->vertexDualArea(he.tipVertex()) * localLapH2 *
        vertexNormal2;
    count++;
  }
}

void System::globalUpdateAfterMutation() {
  // Update the distribution matrix when topology changes
  if (P.eta != 0) {
    D = vpg->d0.transpose().cwiseAbs() / 2;
    // D = vpg->d0.transpose();
    // for (int k = 0; k < D.outerSize(); ++k) {
    //   for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it) {
    //     it.valueRef() = 0.5;
    //   }
    // }
  }

  // Update mask when topology changes
  if (O.isOpenMesh) {
    mask.fill(true);
    boundaryMask(*mesh, mask);
    for (gcs::Vertex v : mesh->vertices()) {
      if (!mask[v]) {
        vpg->inputVertexPositions[v].z = 0;
      }
    }
  }

  // Update spontaneous curvature and bending rigidity when topology changes
  if (!O.isLocalCurvature) {
    H0.raw().setConstant(mesh->nVertices(), 1, P.H0);
    Kb.raw().setConstant(mesh->nVertices(), 1, P.Kb);
  }

  // Update the vertex when topology changes
  if (!O.isFloatVertex) {
    for (gcs::Vertex v : mesh->vertices()) {
      if (thePointTracker[v]) {
        thePoint = gcs::SurfacePoint(v);
      }
    }
    if (thePointTracker.raw().cast<int>().sum() != 1) {
      throw std::runtime_error("globalUpdateAfterMutation: there is no "
                               "unique/existing \"the\" point!");
    }
  }
}

} // namespace mem3dg
