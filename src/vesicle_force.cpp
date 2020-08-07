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
//     Padmini Rangmani (prangamani@eng.ucsd.edu)
//

#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include "ddgsolver/force.h"
#include "ddgsolver/meshops.h"
#include "ddgsolver/util.h"

namespace ddgsolver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getConservativeForces() {
  
  /// 0. GENERAL
  // map the VertexData XXForces to eigen matrix XXForces_e
  auto bendingForces_e = EigenMap<double, 3>(bendingForces);
  auto pressureForces_e = EigenMap<double, 3>(pressureForces);
  auto stretchingForces_e = EigenMap<double, 3>(stretchingForces);

  /// A. BENDING FORCE
  // update the (tufted) mass and conformal Laplacian matrix
  if (isTuftedLaplacian) {
    getTuftedLaplacianAndMass(M, L, mesh, vpg, mollifyFactor);
  } else {
    M = vpg.vertexLumpedMassMatrix;
    L = vpg.cotanLaplacian;
  }
  // Cache the inverse mass matrix
  M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();
  // Gaussian curvature per vertex Area
  auto &KG = vpg.vertexGaussianCurvatures.raw();
  // number of vertices for convenience
  std::size_t n_vertices = (mesh.nVertices());
  // map ivp to eigen matrix position
  auto positions = EigenMap<double, 3>(vpg.inputVertexPositions);
  // the build-in angle-weighted vertex normal
  auto vertexAngleNormal_e = EigenMap<double, 3>(vpg.vertexNormals);
  // calculate mean curvature
  H = rowwiseDotProduct(L * positions / 2.0, vertexAngleNormal_e);
  // calculate the Laplacian of mean curvature H
  Eigen::Matrix<double, Eigen::Dynamic, 1> lap_H = L * M_inv * H;
  // initialize the spontaneous curvature matrix
  H0.setConstant(n_vertices, 1, P.H0);
  // initialize and calculate intermediary result scalarTerms
  Eigen::Matrix<double, Eigen::Dynamic, 1> scalarTerms =
      M_inv * rowwiseProduct(H, H) + rowwiseProduct(H, H0) - KG;
  // initialize and calculate intermediary result productTerms
  Eigen::Matrix<double, Eigen::Dynamic, 1> productTerms;
  productTerms.resize(n_vertices, 3);
  productTerms = rowwiseProduct(scalarTerms, M_inv * H - H0);
  // calculate bendingForce
  bendingForces_e = -2.0 * P.Kb * rowwiseScaling(productTerms + lap_H, vertexAngleNormal_e);

  /// B. PRESSURE FORCES
  volume = 0;
  double face_volume;
  gcs::FaceData<int> sign_of_volume(mesh);
  for (gcs::Face f : mesh.faces()) {
    face_volume = signedVolumeFromFace(f, vpg);
    volume += face_volume;
    if (face_volume < 0) {
      sign_of_volume[f] = -1;
    } else {
      sign_of_volume[f] = 1;
    }
  } 
  pressureForces_e =
      rowwiseScaling(-2.0 * P.Kv * (volume - maxVolume * P.Vt) /
                         (maxVolume * P.Vt) * vpg.vertexDualAreas.raw(),
                     vertexAngleNormal_e); 


  /// C. STRETCHING FORCES (CAPILLARY PRESSURE)
  const gcs::FaceData<gc::Vector3> &face_n = vpg.faceNormals;
  const gcs::FaceData<double> &face_a = vpg.faceAreas;
  Eigen::Matrix<double, Eigen::Dynamic, 1> faceArea_e = vpg.faceAreas.raw();
  surfaceArea = faceArea_e.sum();
  stretchingForces_e =
      rowwiseScaling(-2.0 * P.Ksg * (surfaceArea - targetSurfaceArea) /
                         targetSurfaceArea * 2.0 * H,
                     vertexAngleNormal_e);

  /// D. LOCAL REGULARIZATION
  if ((P.Ksl != 0) && (P.Kse != 0)) {
    for (gcs::Vertex v : mesh.vertices()) {

      for (gcs::Halfedge he : v.outgoingHalfedges()) {
        gcs::Halfedge base_he = he.next();

        // Stretching forces
        gc::Vector3 edgeGradient = -vecFromHalfedge(he, vpg).normalize();
        gc::Vector3 base_vec = vecFromHalfedge(base_he, vpg);
        gc::Vector3 gradient = -gc::cross(base_vec, face_n[he.face()]);
        assert((gc::dot(gradient, vecFromHalfedge(he, vpg))) < 0);
        if (P.Ksl != 0) {
          stretchingForces[v] +=
              -2.0 * P.Ksl * gradient *
              (face_a[base_he.face()] - targetFaceAreas[base_he.face()]) /
              targetFaceAreas[base_he.face()];
        }

        if (P.Kse != 0) {
          stretchingForces[v] +=
              -P.Kse * edgeGradient *
              (vpg.edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
              targetEdgeLengths[he.edge()];
        }
      }
    }
  }

}
} // end namespace ddgsolver
