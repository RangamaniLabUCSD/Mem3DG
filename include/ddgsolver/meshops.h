
#pragma once

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/ply_halfedge_mesh_data.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include "ddgsolver/util.h"
#include "ddgsolver/macros.h"

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
signedVolumeFromFace(gcs::Face &f, gcs::VertexPositionGeometry &vpg) {
  gc::Vector3 p[3];
  size_t i = 0;
  for (gcs::Vertex v : f.adjacentVertices()) {
    p[i] = vpg.inputVertexPositions[v];
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
DLL_PUBLIC inline Eigen::Matrix<double, Eigen::Dynamic, 1> rowwiseDotProduct
(Eigen::Matrix<double, Eigen::Dynamic, 3> A, Eigen::Matrix<double, Eigen::Dynamic, 3> B) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> C =
    ((A.array() * B.array()).rowwise().sum()).matrix();
  return C;
}

/**
 * @brief helper function for rowwise scaling of a matrix using a vector s
 *
 * @param Eigen vector a
 * @param Eigen matrix B
 * @return Eigen matrix C
 */
/// 
DLL_PUBLIC inline Eigen::Matrix<double, Eigen::Dynamic, 3> rowwiseScaling
(Eigen::Matrix<double, Eigen::Dynamic, 1> a, Eigen::Matrix<double, Eigen::Dynamic, 3> B) {
  Eigen::Matrix<double, Eigen::Dynamic, 3> C =
    (B.array().colwise() * a.array()).matrix();
  return C;
}
} // namespace ddgsolver
