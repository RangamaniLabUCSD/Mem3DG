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

#pragma once

#include "geometrycentral/surface/halfedge_factories.h"
#include "geometrycentral/surface/tufted_laplacian.h"
#include <geometrycentral/surface/surface_mesh.h>

#include <Eigen/Core>

#include <math.h>

#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/util.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

/**
 * @brief Get volume from a face
 *
 * @param f
 * @param vpg
 * @return double
 */
DLL_PUBLIC inline double
signedVolumeFromFace(gcs::Face &f, gcs::VertexPositionGeometry &vpg, gc::Vector3 center) {
  gc::Vector3 p[3];
  size_t i = 0;
  for (gcs::Vertex v : f.adjacentVertices()) {
    p[i] = vpg.inputVertexPositions[v] - center;
    i++;
  }
  double v321 = p[2].x * p[1].y * p[0].z;
  double v231 = p[1].x * p[2].y * p[0].z;
  double v312 = p[2].x * p[0].y * p[1].z;
  double v132 = p[0].x * p[2].y * p[1].z;
  double v213 = p[1].x * p[0].y * p[2].z;
  double v123 = p[0].x * p[1].y * p[2].z;

  return (-v321 + v231 + v312 - v132 - v213 + v123) / 6.0;
}

/**
 * @brief Get the vector from halfedge vertices
 *
 * @param he
 * @param vpg
 * @return gc::Vector3
 */
DLL_PUBLIC inline gc::Vector3
vecFromHalfedge(gcs::Halfedge &he, gcs::VertexPositionGeometry &vpg) {
  return vpg.inputVertexPositions[he.next().vertex()] -
         vpg.inputVertexPositions[he.vertex()];
}

/**
 * @brief helper function for taking rowwise dot product of two matrices
 *
 * @param Eigen matrix A
 * @param Eigen matrix B
 * @return Eigen matrix V
 */
DLL_PUBLIC inline Eigen::Matrix<double, Eigen::Dynamic, 1>
rowwiseDotProduct(Eigen::Matrix<double, Eigen::Dynamic, 3> A,
                  Eigen::Matrix<double, Eigen::Dynamic, 3> B) {
  return ((A.array() * B.array()).rowwise().sum()).matrix();
}

/**
 * @brief helper function for taking rowwise cross product of two matrices
 *
 * @param Eigen matrix A
 * @param Eigen matrix B
 * @return Eigen matrix V
 */
DLL_PUBLIC inline Eigen::Matrix<double, Eigen::Dynamic, 3>
rowwiseCrossProduct(Eigen::Matrix<double, Eigen::Dynamic, 3> A,
                  Eigen::Matrix<double, Eigen::Dynamic, 3> B) {
  Eigen::Matrix<double, Eigen::Dynamic, 3> C;
  if (A.rows() != B.rows()) {
    throw std::runtime_error("The input matrices must have same sizes!");
    return C;
  } else {
    C.resize(A.rows(), 3);
    for (size_t i = 0; i < A.rows(); i++) {
      C.row(i) = A.row(i).cross(B.row(i));
    }
    return C;
  }
}

/**
 * @brief helper function for taking rowwise product of two vectors
 *
 * @param Eigen vector A
 * @param Eigen vector B
 * @return Eigen matrix V
 */
DLL_PUBLIC inline Eigen::Matrix<double, Eigen::Dynamic, 1>
rowwiseProduct(Eigen::Matrix<double, Eigen::Dynamic, 1> A,
               Eigen::Matrix<double, Eigen::Dynamic, 1> B) {
  return (A.array() * B.array()).matrix();
}

/**
 * @brief helper function for rowwise scaling of a matrix using a vector s
 *
 * @param Eigen vector a
 * @param Eigen matrix B
 * @return Eigen matrix C
 */
DLL_PUBLIC inline Eigen::Matrix<double, Eigen::Dynamic, 3>
rowwiseScaling(Eigen::Matrix<double, Eigen::Dynamic, 1> a,
               Eigen::Matrix<double, Eigen::Dynamic, 3> B) {
  return (B.array().colwise() * a.array()).matrix();
}

/**
 * @brief helper function for constructing tufted laplacian and mass matrix
 *
 * @param mass matrix M
 * @param Laplacian matrix L
 * @param surfaceMesh mesh
 * @param vertexPositionGeometry vpg
 * @param double mollifyFactor
 */
DLL_PUBLIC inline void getTuftedLaplacianAndMass(
    Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &L,
    gcs::SurfaceMesh &mesh, gcs::VertexPositionGeometry &vpg,
    double mollifyFactor) {

  std::vector<gc::Vector3> vecPosition(
      vpg.inputVertexPositions.raw().data(),
      vpg.inputVertexPositions.raw().data() +
          vpg.inputVertexPositions.raw().size());

  std::unique_ptr<gcs::SurfaceMesh> generalMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> generalVpg;

  std::tie(generalMesh, generalVpg) = gcs::makeGeneralHalfedgeAndGeometry(
      mesh.getFaceVertexList(), vecPosition);
  std::tie(L, M) =
      gcs::buildTuftedLaplacian(*generalMesh, *generalVpg, mollifyFactor);
}

/**
 * @brief Apply vertex shift by moving the vertices chosen for integration to
 * the Barycenter of the it neighbors
 *
 * @param surfaceMesh mesh
 * @param vertexPositionGeometry vpg
 * @param Eigen boolean vector mask
 */
DLL_PUBLIC inline void
vertexShift(gcs::SurfaceMesh &mesh, gcs::VertexPositionGeometry &vpg,
            Eigen::Matrix<bool, Eigen::Dynamic, 1> mask) {

  for (gcs::Vertex v : mesh.vertices()) {

    if (mask(vpg.vertexIndices[v])) {

      gc::Vector3 baryCenter{0.0, 0.0, 0.0};
      double n_vAdj = 0.0;
      for (gcs::Vertex vAdj : v.adjacentVertices()) {
        baryCenter += vpg.inputVertexPositions[vAdj];
        n_vAdj += 1.0;
      }
      baryCenter /= n_vAdj;
      for (gcs::Halfedge he : v.outgoingHalfedges()) {
        gcs::Halfedge base_he = he.next();
        vpg.inputVertexPositions[v] =
            baryCenter - gc::dot(vpg.vertexNormals[v],
                                 baryCenter - vpg.inputVertexPositions[v]) *
                             vpg.vertexNormals[v];
      }
    }
  }
}

/**
 * @brief Remove the rigid body translation
 *
 * @param Eigen pressure matrix
 */
DLL_PUBLIC inline void 
 removeTranslation(Eigen::Matrix<double, Eigen::Dynamic, 3> &pressure) {
  pressure = pressure.rowwise() -
                     ((pressure).colwise().sum() / pressure.rows());
}

/**
 * @brief Remove the rigid body Rotation
 *
 * @param Eigen pressure matrix
 * @param Eigen position matrix
 */
DLL_PUBLIC inline void removeRotation(
    Eigen::Matrix<double, Eigen::Dynamic, 3>
        position, Eigen::Matrix<double, Eigen::Dynamic, 3> &pressure) {
  pressure =
      pressure.rowwise() -
      (rowwiseCrossProduct(position,
                           pressure)
           .colwise()
           .sum() /
       pressure.rows());
}

/**
 * @brief Gaussian distribution
 *
 * @param distance vector 
 * @param standard deviation
 */
DLL_PUBLIC inline void
gaussianDistribution(Eigen::Matrix<double, Eigen::Dynamic, 1>& distribution,
                         Eigen::Matrix<double, Eigen::Dynamic, 1> distance,
                     double stdDev) {
  distribution = (-distance.array() * distance.array() / (2 * stdDev * stdDev)).exp() /
      (stdDev * pow(M_PI * 2, 0.5));
}

/**
 * @brief height = 1 tanh step function with radius r
 *
 * @param (double) sharpness of transition
 * @param (double) radius of height = 1
 * @param (Eigen vector) distance vector
 * 
 */
DLL_PUBLIC inline void
tanhDistribution(Eigen::Matrix<double, Eigen::Dynamic, 1>& distribution,
                     Eigen::Matrix<double, Eigen::Dynamic, 1> distance,
                     double sharpness, double radius) {
  distribution.resize(distance.rows(), 1);
  for (size_t i = 0; i < distance.rows(); i++) {
    distribution[i] = 0.5 * (1 + tanh(sharpness * (radius - distance[i])));
  }
}

/**
 * @brief Get the Face Vertex Matrix object
 *
 * @param mesh   Mesh of interest
 * @return Eigen matrix of uint32_t indices
 */
DLL_PUBLIC inline Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
getFaceVertexMatrix(gcs::SurfaceMesh &mesh) {
  if (!mesh.isTriangular())
    throw std::runtime_error("Mesh is not triangular.");

  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> result(mesh.nFaces(), 3);

  gcs::VertexData<std::size_t> vInd = mesh.getVertexIndices();
  
  std::size_t i = 0;
  for (gcs::Face f : mesh.faces()) {
    std::uint32_t a, b, c;
    gcs::Halfedge tmp = f.halfedge();
    a = vInd[tmp.vertex()];
    tmp = tmp.next();
    b = vInd[tmp.vertex()];
    tmp = tmp.next();
    c = vInd[tmp.vertex()];

    result.row(i++) << a, b, c;
  }
  return result;
}

} // namespace ddgsolver
