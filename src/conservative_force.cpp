
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
    /// A. BENDING FORCE
    // Initialize the mass matrix
    M = vpg.vertexLumpedMassMatrix;
    M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();
    // Initialize the conformal Laplacian matrix
    L = vpg.cotanLaplacian;

    // Gaussian curvature per vertex Area
    Eigen::Matrix<double, Eigen::Dynamic, 1> KG =
      M_inv * (vpg.vertexGaussianCurvatures.toMappedVector());

    // number of vertices for convenience
    std::size_t n_vertices = (mesh.nVertices());

    // map ivp to eigen matrix position
    auto positions = ddgsolver::EigenMap<double, 3>(vpg.inputVertexPositions);

    // map the VertexData bendingForces to eigen matrix bendingForces_e
    auto bendingForces_e = ddgsolver::EigenMap<double, 3>(bendingForces);
    bendingForces_e.setZero();

    // the build-in angle-weighted vertex normal
    auto vertexAngleNormal_e = ddgsolver::EigenMap<double, 3>(vpg.vertexNormals);

    // calculate mean curvature per vertex area and map it to angle-weighted normal
    Hn = rowwiseScaling(rowwiseDotProduct(vertexAngleNormal_e,
      M_inv * L * positions / 2.0), vertexAngleNormal_e);

    // calculate the Laplacian of mean curvature H 
    Eigen::Matrix<double, Eigen::Dynamic, 3> lap_H = M_inv * L * Hn;

    // initialize the spontaneous curvature matrix
    H0n = H0 * vertexAngleNormal_e;

    // initialize and calculate intermediary result scalerTerms, set to zero if negative
    Eigen::Matrix<double, Eigen::Dynamic, 1> scalerTerms =
      rowwiseDotProduct(Hn, Hn) + rowwiseDotProduct(H0n, H0n) - KG;
    Eigen::Matrix<double, Eigen::Dynamic, 1> zeroMatrix;
    zeroMatrix.resize(n_vertices, 1);
    zeroMatrix.setZero();
    scalerTerms = scalerTerms.array().max(zeroMatrix.array());

    // initialize and calculate intermediary result productTerms
    Eigen::Matrix<double, Eigen::Dynamic, 3> productTerms;
    productTerms.resize(n_vertices, 3);
    productTerms = 2 * rowwiseScaling(scalerTerms, Hn - H0n);

    // calculate bendingForce
    bendingForces_e = M * (-2.0 * Kb * (productTerms + lap_H));

    /// B. PRESSURE FORCES
    pressureForces.fill({ 0.0, 0.0, 0.0 });
    volume = 0;
    double face_volume;
    gcs::FaceData<int> sign_of_volume(mesh);
    for (gcs::Face f : mesh.faces()) {
      face_volume = signedVolumeFromFace(f, vpg);
      volume += face_volume;
      if (face_volume < 0) {
        sign_of_volume[f] = -1;
      }
      else {
        sign_of_volume[f] = 1;
      }
    }

    /// C. STRETCHING FORCES
    stretchingForces.fill({ 0.0,0.0,0.0 });
    const gcs::FaceData<gc::Vector3>& face_n = vpg.faceNormals;
    const gcs::FaceData<double>& face_a = vpg.faceAreas;
    auto faceArea_e = EigenMap(vpg.faceAreas);
    surfaceArea = faceArea_e.sum();

    /// D. LOOPING VERTICES
    for (gcs::Vertex v : mesh.vertices()) {

      for (gcs::Halfedge he : v.outgoingHalfedges()) {

        // Pressure forces
        gcs::Halfedge base_he = he.next();
        gc::Vector3 p1 = vpg.inputVertexPositions[base_he.vertex()];
        gc::Vector3 p2 = vpg.inputVertexPositions[base_he.next().vertex()];
        gc::Vector3 dVdx = 0.5 * gc::cross(p1, p2) / 3.0;
        assert(gc::dot(dVdx, vpg.inputVertexPositions[v] - p1) *
          sign_of_volume[he.face()] >
          0);
        pressureForces[v] +=
          -0.5 * Kv * (volume - maxVolume * Vt) / (maxVolume * Vt) * dVdx;

        // Stretching forces
        gc::Vector3 edgeGradient = -vecFromHalfedge(he, vpg).normalize();
        gc::Vector3 base_vec = vecFromHalfedge(base_he, vpg);
        gc::Vector3 gradient = -gc::cross(base_vec, face_n[he.face()]);
        assert((gc::dot(gradient, vecFromHalfedge(he, vpg))) < 0);
        if (Ksl != 0) {
          stretchingForces[v] += -2 * Ksl * gradient *
            (face_a[base_he.face()] - initialFaceAreas[base_he.face()]) /
            initialFaceAreas[base_he.face()];
        }
        if (Ksg != 0) {
          stretchingForces[v] +=
            -2 * Ksg * gradient * (surfaceArea - initialSurfaceArea) / initialSurfaceArea;
        }
        if (Kse != 0) {
          stretchingForces[v] += -Kse * edgeGradient *
            (vpg.edgeLengths[he.edge()] - targetEdgeLength[he.edge()]) / targetEdgeLength[he.edge()];
        }
      
      }
    }

  }
} // end namespace ddgsolver
