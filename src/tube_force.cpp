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
#include <geometrycentral/utilities/vector3.h>
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"

#include <Eigen/Core>

#include "mem3dg/solver/force.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/util.h"

namespace ddgsolver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getTubeForces() {

  /// A. BENDING FORCE
  if (P.Kb != 0) {
    
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

    // map the VertexData bendingForces to eigen matrix bendingForces_e
    auto bendingForces_e = EigenMap<double, 3>(bendingForces);

    // the build-in angle-weighted vertex normal
    auto vertexAngleNormal_e = EigenMap<double, 3>(vpg.vertexNormals);

    // calculate mean curvature
    H = rowwiseDotProduct(L * positions / 2.0, vertexAngleNormal_e);

    // calculate the Laplacian of mean curvature H
    Eigen::Matrix<double, Eigen::Dynamic, 3> lap_H =
        L * M_inv * rowwiseScaling(H, vertexAngleNormal_e);

    // initialize the spontaneous curvature matrix
    H0.setConstant(n_vertices, 1, P.H0);

    // initialize and calculate intermediary result scalarTerms, set to zero if
    // negative
    Eigen::Matrix<double, Eigen::Dynamic, 1> scalarTerms =
        M_inv * rowwiseProduct(H, H) + rowwiseProduct(H, H0) - KG;
    /*Eigen::Matrix<double, Eigen::Dynamic, 1> zeroMatrix;
    zeroMatrix.resize(n_vertices, 1);
    zeroMatrix.setZero();
    scalarTerms = scalarTerms.array().max(zeroMatrix.array());*/

    // initialize and calculate intermediary result productTerms
    Eigen::Matrix<double, Eigen::Dynamic, 3> productTerms;
    productTerms.resize(n_vertices, 3);
    productTerms =
        2 * rowwiseScaling(rowwiseProduct(scalarTerms, M_inv * H - H0),
                           vertexAngleNormal_e);

    // calculate bendingForce
    bendingForces_e = -2.0 * P.Kb * (productTerms + lap_H);
  }

  /// B. PRESSURE FORCES
  pressureForces.fill({0.0, 0.0, 0.0});
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

  /// C. STRETCHING FORCES
  stretchingForces.fill({0.0, 0.0, 0.0});
  const gcs::FaceData<gc::Vector3> &face_n = vpg.faceNormals;
  const gcs::FaceData<double> &face_a = vpg.faceAreas;
  Eigen::Matrix<double, Eigen::Dynamic, 1> faceArea_e = vpg.faceAreas.raw();
  surfaceArea = faceArea_e.sum();

  /// D. LOOPING VERTICES
  for (gcs::Vertex v : mesh.vertices()) {

    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gcs::Halfedge base_he = he.next();

      // Pressure forces = 0 for patch tube simulation
      if (P.Kv != 0) {

        P.Kv = 0;
        std::cout << "\n"
                  << "Kv suppose to be 0 in patch tube simulation!!!!"
                  << "\n"
                  << "Kv = 0 now" << std::endl;

        gc::Vector3 p1 = vpg.inputVertexPositions[base_he.vertex()];
        gc::Vector3 p2 = vpg.inputVertexPositions[base_he.next().vertex()];
        gc::Vector3 dVdx = 0.5 * gc::cross(p1, p2) / 3.0;
        assert(gc::dot(dVdx, vpg.inputVertexPositions[v] - p1) *
                   sign_of_volume[he.face()] >
               0);
        pressureForces[v] += -2.0 * P.Kv * (volume - maxVolume * P.Vt) /
                             (maxVolume * P.Vt) * dVdx;
      }

      // Stretching forces
      gc::Vector3 edgeGradient = -vecFromHalfedge(he, vpg).normalize();
      gc::Vector3 base_vec = vecFromHalfedge(base_he, vpg);
      gc::Vector3 gradient = -gc::cross(base_vec, face_n[he.face()]);
      assert((gc::dot(gradient, vecFromHalfedge(he, vpg))) < 0);
      
      // patch simulation assumes constant surface tension
      if (P.Ksl != 0) {
        stretchingForces[v] +=
            -2.0 * P.Ksl * gradient;
      }
      if (P.Ksg != 0) {
        stretchingForces[v] += -2.0 * P.Ksg * gradient;
      }

      // the cubic penalty is for regularizing the mesh,
      // need better physical interpretation or alternative method
      if (P.Kse != 0) {
        double strain =
            (vpg.edgeLengths[he.edge()] - targetEdgeLengths[he.edge()]) /
            targetEdgeLengths[he.edge()];
        stretchingForces[v] +=
            -P.Kse * edgeGradient * strain * strain * strain;
      }
    }
  }
}
} // end namespace ddgsolver
