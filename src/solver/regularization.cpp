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
#include "mem3dg/constants.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/system.h"
#include <Eigen/Core>
#include <cmath>

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::computeRegularizationForce() {
  // Note in regularization, it is preferred to use immediate calculation rather
  // than cached one
  for (gcs::Vertex v : mesh->vertices()) {
    if (!v.isBoundary()) {
      for (gcs::Halfedge he : v.outgoingHalfedges()) {
        gcs::Edge e = he.edge();
        // Conformal regularization
        if (meshProcessor.meshRegularizer.Kst != 0 && !e.isBoundary()) {
          gcs::Halfedge jl = he.next();
          gcs::Halfedge li = jl.next();
          gcs::Halfedge ik = he.twin().next();
          gcs::Halfedge kj = ik.next();

          gc::Vector3 grad_li = vecFromHalfedge(li, *vpg).normalize();
          gc::Vector3 grad_ik = vecFromHalfedge(ik.twin(), *vpg).normalize();
          forces.regularizationForce[v] +=
              -meshProcessor.meshRegularizer.Kst *
              (meshProcessor.meshRegularizer.computeLengthCrossRatio(*vpg, e) -
               meshProcessor.meshRegularizer.refLcrs[e.getIndex()]) /
              meshProcessor.meshRegularizer.refLcrs[e.getIndex()] *
              (vpg->edgeLength(kj.edge()) / vpg->edgeLength(jl.edge())) *
              (grad_li * vpg->edgeLength(ik.edge()) -
               grad_ik * vpg->edgeLength(li.edge())) /
              vpg->edgeLength(ik.edge()) / vpg->edgeLength(ik.edge());
        }

        // Local area regularization
        if (meshProcessor.meshRegularizer.Ksl != 0 && he.isInterior()) {
          gcs::Halfedge base_he = he.next();
          gc::Vector3 base_vec = vecFromHalfedge(base_he, *vpg);
          gc::Vector3 localAreaGradient =
              -gc::cross(base_vec, vpg->faceNormal(he.face()));
          auto &referenceArea =
              (v.isBoundary()
                   ? meshProcessor.meshRegularizer
                         .refFaceAreas[base_he.face().getIndex()]
                   : meshProcessor.meshRegularizer.meanTargetFaceArea);
          forces.regularizationForce[v] +=
              -meshProcessor.meshRegularizer.Ksl * localAreaGradient *
              (vpg->faceArea(base_he.face()) - referenceArea);
        }

        // local edge regularization
        if (meshProcessor.meshRegularizer.Kse != 0) {
          gc::Vector3 edgeGradient = -vecFromHalfedge(he, *vpg).normalize();
          auto &referenceLength =
              (v.isBoundary()
                   ? meshProcessor.meshRegularizer.refEdgeLengths[e.getIndex()]
                   : meshProcessor.meshRegularizer.meanTargetEdgeLength);
          forces.regularizationForce[v] +=
              -meshProcessor.meshRegularizer.Kse * edgeGradient *
              (vpg->edgeLength(e) - referenceLength);
        }
      }
    }
  }

  // post processing regularization force
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg->vertexNormals);
  auto regularizationForce_e =
      gc::EigenMap<double, 3>(forces.regularizationForce);

  // remove the normal component
  regularizationForce_e -= rowwiseScalarProduct(
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
  //     rowwiseScalarProduct(mask.raw().cast<double>(), regularizationForce_e);
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
    if (gc::sum(forces.forceMask[v]) > 0.5) {
      if (v.isBoundary()) {
        gcs::Vertex v1 = v;
        gcs::Vertex v2 = v;
        gc::Vector3 baryCenter{0.0, 0.0, 0.0};
        int n_vAdj = 0;
        for (gcs::Vertex vAdj : v.adjacentVertices()) {
          if (vAdj.isBoundary()) {
            // std::cout << "v: " << v.getIndex() << std::endl;
            // std::cout << "v1:  " << v1.getIndex() << std::endl;
            // std::cout << "v2: " << v2.getIndex() << std::endl;
            if (v1 == v) {
              v1 = vAdj;
            } else if (v2 == v) {
              v2 = vAdj;
            }
            // v1 = (v1 == v) ? vAdj : v1;
            // v2 = (v2 == v) ? vAdj : v2;
            n_vAdj += 1;
          }
        }
        if (n_vAdj != 2) {
          mem3dg_runtime_error(
              "Number of neighbor vertices on boundary is not 2!");
        }
        baryCenter =
            (vpg->inputVertexPositions[v1] + vpg->inputVertexPositions[v2]) / 2;
        gc::Vector3 faceNormal = gc::cross(
            vpg->inputVertexPositions[v1] - vpg->inputVertexPositions[v],
            vpg->inputVertexPositions[v2] - vpg->inputVertexPositions[v]);
        gc::Vector3 sideNormal =
            gc::cross(faceNormal, vpg->inputVertexPositions[v1] -
                                      vpg->inputVertexPositions[v2])
                .normalize();
        vpg->inputVertexPositions[v] =
            baryCenter -
            gc::dot(sideNormal, baryCenter - vpg->inputVertexPositions[v]) *
                sideNormal;
      } else {
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
}

bool System::edgeFlip() {
  // Note in regularization, it is preferred to use immediate calculation rather
  // than cached one
  bool isFlipped = false;
  gcs::EdgeData<bool> isOrigEdge(*mesh, true);
  // flip edge if not delauney
  for (gcs::Edge e : mesh->edges()) {
    if (!isOrigEdge[e] || e.isBoundary()) {
      continue;
    }
    gcs::Halfedge he = e.halfedge();
    if (gc::sum(forces.forceMask[he.vertex()] +
                forces.forceMask[he.twin().vertex()]) < 0.5) {
      continue;
    }

    if (meshProcessor.meshMutator.ifFlip(e, *vpg)) {
      bool sucess = mesh->flip(e);
      isOrigEdge[e] = false;
      isFlipped = true;
      meshProcessor.meshMutator.markVertices(mutationMarker, he.tailVertex());
      meshProcessor.meshMutator.markVertices(mutationMarker, he.tipVertex());
    }
  }

  if (isFlipped)
    mesh->compress();

  return isFlipped;
}

bool System::growMesh() {
  // Note in regularization, it is preferred to use immediate calculation rather
  // than cached one
  bool isGrown = false;
  int count = 0;
  gcs::EdgeData<bool> isOrigEdge(*mesh, true);
  // gcs::VertexData<bool> isOrigVertex(*mesh, true);

  // expand the mesh when area is too large
  for (gcs::Edge e : mesh->edges()) {

    // don't keep processing new edges
    if (!isOrigEdge[e])
      continue;

    // alias the halfedge
    gcs::Halfedge he = e.halfedge();

    // gather both vertices and their properties
    gcs::Vertex vertex1 = he.tipVertex(), vertex2 = he.tailVertex();
    gc::Vector3 vertex1Pos = vpg->vertexPositions[vertex1];
    gc::Vector3 vertex2Pos = vpg->vertexPositions[vertex2];
    gc::Vector3 vertex1Vel = velocity[vertex1];
    gc::Vector3 vertex2Vel = velocity[vertex2];
    double vertex1GeoDist = geodesicDistanceFromPtInd[vertex1];
    double vertex2GeoDist = geodesicDistanceFromPtInd[vertex2];
    double vertex1Phi = proteinDensity[vertex1];
    double vertex2Phi = proteinDensity[vertex2];
    gc::Vector3 vertex1ForceMask = forces.forceMask[vertex1];
    gc::Vector3 vertex2ForceMask = forces.forceMask[vertex2];
    bool vertex1PointTracker = thePointTracker[vertex1];
    bool vertex2PointTracker = thePointTracker[vertex2];

    // don't keep processing static vertices
    if (gc::sum(vertex1ForceMask + vertex2ForceMask) < 0.5)
      continue;

    // Spltting
    if (meshProcessor.meshMutator.ifSplit(e, *vpg)) {
      count++;
      // split the edge
      gcs::Vertex newVertex = mesh->splitEdgeTriangular(e).vertex();

      // update quantities
      // Note: think about conservation of energy, momentum and angular
      // momentum
      // averageData(vpg->inputVertexPositions, vertex1, vertex2, newVertex);
      // averageData(velocity, vertex1, vertex2, newVertex);
      // averageData(geodesicDistanceFromPtInd, vertex1, vertex2, newVertex);
      // averageData(proteinDensity, vertex1, vertex2, newVertex);
      vpg->vertexPositions[newVertex] = 0.5 * (vertex1Pos + vertex2Pos);
      velocity[newVertex] = 0.5 * (vertex1Vel + vertex2Vel);
      geodesicDistanceFromPtInd[newVertex] =
          0.5 * (vertex1GeoDist + vertex2GeoDist);
      proteinDensity[newVertex] = 0.5 * (vertex1Phi + vertex2Phi);
      thePointTracker[newVertex] = false;
      forces.forceMask[newVertex] = gc::Vector3{1, 1, 1};

      // isOrigVertex[newVertex] = false;
      for (gcs::Edge e : newVertex.adjacentEdges()) {
        isOrigEdge[e] = false;
      }

      meshProcessor.meshMutator.markVertices(mutationMarker, newVertex);
      // mutationMarker[newVertex] = true;

      isGrown = true;
    } else if (meshProcessor.meshMutator.ifCollapse(e, *vpg)) { // Collapsing
      // collapse the edge
      gcs::Vertex newVertex = mesh->collapseEdgeTriangular(e);

      if (newVertex != gcs::Vertex()) {
        count++;
        // update quantities
        // Note: think about conservation of energy, momentum and angular
        // momentum
        vpg->vertexPositions[newVertex] =
            gc::sum(vertex1ForceMask) < 2.5   ? vertex1Pos
            : gc::sum(vertex2ForceMask) < 2.5 ? vertex2Pos
                                              : (vertex1Pos + vertex2Pos) / 2;
        // averageData(velocity, vertex1, vertex2, newVertex);
        // averageData(geodesicDistanceFromPtInd, vertex1, vertex2, newVertex);
        // averageData(proteinDensity, vertex1, vertex2, newVertex);
        velocity[newVertex] = 0.5 * (vertex1Vel + vertex2Vel);
        geodesicDistanceFromPtInd[newVertex] =
            0.5 * (vertex1GeoDist + vertex2GeoDist);
        proteinDensity[newVertex] = 0.5 * (vertex1Phi + vertex2Phi);
        thePointTracker[newVertex] = vertex1PointTracker || vertex2PointTracker;

        // isOrigVertex[newVertex] = false;
        for (gcs::Edge e : newVertex.adjacentEdges()) {
          isOrigEdge[e] = false;
        }

        meshProcessor.meshMutator.markVertices(mutationMarker, newVertex);

        isGrown = true;
      }
    }
  }
  if (isGrown)
    mesh->compress();
  return isGrown;
}

void System::mutateMesh(size_t nRepetition) {
  for (size_t i = 0; i < nRepetition; ++i) {
    bool isGrown = false, isFlipped = false;
    mutationMarker.fill(false);

    // vertex shift for regularization
    if (meshProcessor.meshMutator.shiftVertex) {
      vertexShift();
    }

    // split edge and collapse edge
    if (meshProcessor.meshMutator.isSplitEdge ||
        meshProcessor.meshMutator.isCollapseEdge) {
      isGrown = isGrown || growMesh();
    }

    // linear edge flip for non-Delauney triangles
    if (meshProcessor.meshMutator.isEdgeFlip) {
      isFlipped = edgeFlip();
      isFlipped = edgeFlip() || isFlipped;
      isFlipped = edgeFlip() || isFlipped;
    }

    // globally update quantities
    if (isGrown || isFlipped) {
      globalUpdateAfterMutation();
    }
  }
}

Eigen::Matrix<bool, Eigen::Dynamic, 1>
System::smoothenMesh(double initStep, double target, size_t maxIteration) {
  // require nonzero bending rigidity in parameters
  if (Kb.raw().sum() == 0) {
    mem3dg_runtime_error(
        "Bending rigidity has to be nonzero to smoothen mesh!");
  }
  // initialize variables
  double stepSize = initStep;
  double pastGradNorm = 1e10;
  size_t num_iter = 0;
  // compute bending forces
  vpg->refreshQuantities();
  computeMechanicalForces();
  EigenVectorX3dr pastForceVec = toMatrix(forces.bendingForceVec);
  // initialize smoothingMask
  Eigen::Matrix<bool, Eigen::Dynamic, 1> smoothingMask =
      outlierMask(forces.bendingForce.raw(), 0.5);
  isSmooth = (smoothingMask.cast<int>().sum() == 0);
  // initialize gradient and compute exit tolerance
  double gradNorm = computeNorm(toMatrix(forces.bendingForceVec));
  double tol = gradNorm * target;

  while (gradNorm > tol && !isSmooth) {
    if (stepSize < 1e-8 * initStep) {
      mem3dg_runtime_message("smoothing operation diverges!");
      break;
    }
    if (num_iter == maxIteration) {
      mem3dg_runtime_message("smoothing operation exceeds max iteration!");
      break;
    }

    // compute bending force if smoothingMask is true
    vpg->refreshQuantities();
    forces.bendingForceVec.fill({0, 0, 0});
    forces.bendingForce.raw().setZero();
    for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
      if (smoothingMask[i]) {
        computeMechanicalForces(i);
      }
    }
    // compute norm of the bending force
    gradNorm = computeNorm(toMatrix(forces.bendingForceVec));
    // recover the position and cut the step size in half
    if (gradNorm > pastGradNorm) {
      toMatrix(vpg->inputVertexPositions) -= pastForceVec * stepSize;
      stepSize /= 2;
      continue;
    }
    // smoothing step
    vpg->inputVertexPositions += forces.bendingForceVec * stepSize;
    pastGradNorm = gradNorm;
    pastForceVec = toMatrix(forces.bendingForceVec);
    num_iter++;
    
  };

  return smoothingMask;
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

    gcs::Vertex v = he.tailVertex();
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
  // update the velocity
  velocity = forces.maskForce(velocity); // important: velocity interpolation
                                         // contaminate the zero velocity
  if (computeKineticEnergy() != 0) {
    double oldKE = energy.kineticEnergy;
    velocity *= pow(oldKE / computeKineticEnergy(), 0.5);
  }

  // Update mask when topology changes (likely not necessary, just for safety)
  if (isOpenMesh) {
    forces.forceMask.fill({1, 1, 1});
    boundaryForceMask(*mesh, forces.forceMask,
                      parameters.boundary.shapeBoundaryCondition);
    forces.proteinMask.fill(1);
    boundaryProteinMask(*mesh, forces.proteinMask,
                        parameters.boundary.proteinBoundaryCondition);
    // for (gcs::Vertex v : mesh->vertices()) {
    //   if (!mask[v]) {
    //     vpg->inputVertexPositions[v].z = 0;
    //   }
    // }
  }

  // Update the vertex when topology changes
  if (!parameters.point.isFloatVertex) {
    for (gcs::Vertex v : mesh->vertices()) {
      if (thePointTracker[v]) {
        thePoint = gcs::SurfacePoint(v);
      }
    }
    if (thePointTracker.raw().cast<int>().sum() != 1) {
      mem3dg_runtime_error("globalUpdateAfterMutation: there is no "
                           "unique/existing \"the\" point!");
    }
  }

  // // Update the distribution matrix when topology changes
  // if (P.eta != 0) {
  //   D = vpg->d0.transpose().cwiseAbs() / 2;
  //   // D = vpg->d0.transpose();
  //   // for (int k = 0; k < D.outerSize(); ++k) {
  //   //   for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it)
  //   {
  //   //     it.valueRef() = 0.5;
  //   //   }
  //   // }
  // }

  // Update spontaneous curvature and bending rigidity when topology changes
  // if (!O.isHeterogeneous) {
  //   proteinDensity.raw().se
  //   // H0.raw().setConstant(mesh->nVertices(), 1, P.H0);
  //   // Kb.raw().setConstant(mesh->nVertices(), 1, P.Kb);
  // }
}

} // namespace solver
} // namespace mem3dg
