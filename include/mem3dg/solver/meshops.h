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

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/halfedge_factories.h"
#include "geometrycentral/surface/tufted_laplacian.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"
#include "geometrycentral/utilities/vector3.h"
#include <geometrycentral/surface/surface_mesh.h>

#include <Eigen/Core>

#include <math.h>
#ifdef __linux__
#include <sys/time.h>
#endif

#include "mem3dg/solver/constants.h"
#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/util.h"

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

/**
 * @brief Sort Eigen vector from largest to smallest
 * @param vec vector to be sorted
 * @param sorted_vec sorted vector
 * @param ind indices of sorted vector based on original vector
 * @return
 */
DLL_PUBLIC inline void sortVector(const Eigen::VectorXd &vec,
                                  Eigen::VectorXd &sorted_vec,
                                  Eigen::VectorXi &ind) {
  ind = Eigen::VectorXi::LinSpaced(vec.size(), 0,
                                   vec.size() - 1); //[0 1 2 3 ... N-1]
  auto rule = [vec](int i, int j) -> bool { return vec(i) > vec(j); };
  // regular expression, as a predicate of sort
  std::sort(ind.data(), ind.data() + ind.size(), rule);
  // The data member function returns a pointer to
  // the first element of VectorXd, similar to
  // begin()
  sorted_vec.resize(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    sorted_vec(i) = vec(ind(i));
  }
}

/**
 * @brief Sort Eigen vector from largest to smallest (overloaded)
 * @param vec vector to be sorted
 * @param sorted_vec sorted vector
 * @return
 */
DLL_PUBLIC inline void sortVector(const Eigen::VectorXd &vec,
                                  Eigen::VectorXd &sorted_vec) {
  Eigen::VectorXi ind =
      Eigen::VectorXi::LinSpaced(vec.size(), 0,
                                 vec.size() - 1); //[0 1 2 3 ... N-1]
  auto rule = [vec](int i, int j) -> bool { return vec(i) > vec(j); };
  // regular expression, as a predicate of sort
  std::sort(ind.data(), ind.data() + ind.size(), rule);
  // The data member function returns a pointer to
  // the first element of VectorXd, similar to
  // begin()
  sorted_vec.resize(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    sorted_vec(i) = vec(ind(i));
  }
}

/**
 * @brief find the middle index between two bounding indices
 * @param l left bound
 * @param r right bound
 * @return middle index
 */
DLL_PUBLIC inline int findMedianIndex(int l, int r) {
  int n = r - l + 1;
  n = (n + 1) / 2 - 1;
  return n + l;
}

/**
 * @brief find the range of data based on percentile (quatile)
 * @param a raw buffer of vector
 * @param n size of vector
 * @param r upper bound of range
 * @param l lower bound of range
 * @return
 */
DLL_PUBLIC inline void findRange(double *a, int n, double &r, double &l) {
  // Index of median of entire data
  int mid_index = findMedianIndex(0, n);

  // first 1/128
  r = a[findMedianIndex(
      0, findMedianIndex(
             0, findMedianIndex(
                    0, findMedianIndex(
                           0, findMedianIndex(
                                  0, findMedianIndex(0, mid_index))))))];

  // last 1/128
  l = a[findMedianIndex(
      findMedianIndex(
          findMedianIndex(
              findMedianIndex(
                  findMedianIndex(findMedianIndex(mid_index + 1, n), n), n),
              n),
          n),
      n)];
}

/**
 * @brief test whether exist outliers to the set of data based range function
 * @param vec data vector
 * @param threshold coefficient used to bound the outlier.
 * For example: (Outlier <--> r) < [threshold * (l <--> r)]
 * @return
 */
DLL_PUBLIC inline bool hasOutlier(const Eigen::VectorXd &vec,
                                  double threshold = 0.5) {
  Eigen::VectorXd sorted_vec;
  sortVector(vec, sorted_vec);
  double r, l, range;
  findRange(sorted_vec.data(), sorted_vec.size(), r, l);
  range = r - l;
  return (sorted_vec[0] - r > threshold * range) ||
         (l - sorted_vec[sorted_vec.size() - 1] > threshold * range);
}

/**
 * @brief Signal handler for pybindf
 */
DLL_PUBLIC inline void signalHandler(int signum) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";
  exit(signum);
}

// /**
//  * @brief close a open mesh
//  *
//  * @param f
//  * @param vpg
//  * @return double
//  */
// DLL_PUBLIC inline std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
//                              std::unique_ptr<gcs::VertexPositionGeometry>>
// toClosedMesh(
//     Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> &&coords,
//     Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> &&faces,
//     gcs::BoundaryLoop &bl, gc::Vector3 center) {

//   Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> newFaces;
//   newFaces.resize(bl.degree() - 2 + faces.rows());

//   std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
//   std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;

//   std::size_t count = 0;
//   gcs::Vertex theVertex = bl.halfedge().tailVertex();
//   for (gcs::Halfedge e : bl.adjacentHalfedges()) {
//     if (e.tailVertex() != theVertex && e.tipVertex() != theVertex) {
//       newFaces.row(count) << e.tailVertex().getIndex(),
//           e.tipVertex().getIndex(), theVertex.getIndex();
//       count++;
//     }
//   }

//   newFaces.tail(faces.rows() * 3) << faces;

//   return gcs::makeManifoldSurfaceMeshAndGeometry(coords, newFaces);
// }

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
 * @brief Get volume from a face
 *
 * @param f
 * @param vpg
 * @return double
 */
DLL_PUBLIC inline double
signedVolumeFromFace(std::vector<std::size_t> f,
                     gcs::VertexPositionGeometry &vpg) {
  gc::Vector3 p[3];
  size_t i = 0;
  for (std::size_t v : f) {
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
 * @brief Get mesh volume
 *
 * @param f
 * @param vpg
 * @return double
 */
DLL_PUBLIC inline double getMeshVolume(gcs::ManifoldSurfaceMesh &mesh,
                                       gcs::VertexPositionGeometry &vpg,
                                       bool isFillHole = true) {
  double volume = 0;
  for (gcs::Face f : mesh.faces()) {
    volume += signedVolumeFromFace(f, vpg);
  }

  // fill hole for open mesh
  if (mesh.hasBoundary()) {
    if (isFillHole) {
      for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
        gcs::Vertex theVertex = bl.halfedge().tailVertex();
        for (gcs::Halfedge e : bl.adjacentHalfedges()) {
          if (e.tailVertex() != theVertex && e.tipVertex() != theVertex) {
            volume += signedVolumeFromFace(
                std::vector<std::size_t>{e.tailVertex().getIndex(),
                                         e.tipVertex().getIndex(),
                                         theVertex.getIndex()},
                vpg);
          }
        }
      }
    } else {
      throw std::runtime_error("getMeshVolume: mesh is opened, not able to "
                               "compute enclosed volume unless filled holes!");
    }
  }
  return volume;
}

/**
 * @brief Get average data
 *
 * @param meshData
 * @param element1
 * @param element2
 * @param newElement
 * @return
 */
template <typename E, typename T>
DLL_PUBLIC inline void averageData(gc::MeshData<E, T> &meshData,
                                   const E &element1, const E &element2,
                                   const E &newElement) {
  meshData[newElement] = (meshData[element1] + meshData[element2]) / 2;
}

/**
 * @brief get barycentric coordinate from cartesian coordinate
 *
 * @param meshData
 * @param element1
 * @param element2
 * @param newElement
 * @return
 */
DLL_PUBLIC inline gc::Vector3 cartesianToBarycentric(gc::Vector2 &v1,
                                                     gc::Vector2 &v2,
                                                     gc::Vector2 &v3,
                                                     gc::Vector2 &v) {
  gc::Vector3 lambda;

  lambda.x = ((v2.y - v3.y) * (v.x - v3.x) + (v3.x - v2.x) * (v.y - v3.y)) /
             ((v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y));
  lambda.y = ((v3.y - v1.y) * (v.x - v3.x) + (v1.x - v3.x) * (v.y - v3.y)) /
             ((v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y));
  lambda.z = 1 - lambda.x - lambda.y;

  return lambda;
}

/**
 * @brief get coorresponding barycentric coordinate to a face
 *
 * @param baryCoords reference to Barycentric coordinate
 * @param firstHalfedge reference to the halfedge associated with the first
 * vertex of the Barycentric coordinate
 * @return
 */
DLL_PUBLIC inline gc::Vector3
correspondBarycentricCoordinates(gc::Vector3 &baryCoords_,
                                 gcs::Halfedge &firstHalfedge) {
  size_t vertexInd = 0;
  gc::Vector3 baryCoords;
  for (gcs::Vertex v : firstHalfedge.face().adjacentVertices()) {
    if (v == firstHalfedge.vertex()) {
      baryCoords[vertexInd] = baryCoords_.x;
    } else if (v == firstHalfedge.next().vertex()) {
      baryCoords[vertexInd] = baryCoords_.y;
    } else if (v == firstHalfedge.next().next().vertex()) {
      baryCoords[vertexInd] = baryCoords_.z;
    } else if (firstHalfedge.face().isBoundaryLoop()) {
      throw std::runtime_error(
          "correspondBarycentricCoordinates: face is on boundary!");
    } else {
      throw std::runtime_error(
          "correspondBarycentricCoordinates: undefined behavior!");
    }
    vertexInd++;
  }
  return baryCoords;
}

/**
 * @brief Get the vector from halfedge vertices
 *
 * @param he
 * @param vpg
 * @return gc::Vector3
 */
DLL_PUBLIC inline gc::Vector3
vecFromHalfedge(const gcs::Halfedge &he,
                const gcs::VertexPositionGeometry &vpg) {
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
 * @brief helper function for computing the polygon area enclosed by a boundary
 * loop on a mesh
 * @param bl boundary loop on a mesh
 * @param inputVertexPosition embedded vertex position of the mesh
 * @return enclosed polygon area
 */
DLL_PUBLIC inline double
computePolygonArea(const gcs::BoundaryLoop &bl,
                   const gcs::VertexData<gc::Vector3> &inputVertexPositions) {
  gc::Vector3 signedArea{0, 0, 0};
  for (gcs::Halfedge he : bl.adjacentHalfedges()) {
    signedArea += cross(inputVertexPositions[he.tailVertex()],
                        inputVertexPositions[he.tipVertex()]);
  }

  // could be used to poject onto other direction if needed
  return 0.5 * signedArea.norm();
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
 */
DLL_PUBLIC inline void getTuftedLaplacianAndMass(
    Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &L,
    gcs::SurfaceMesh &mesh, gcs::VertexPositionGeometry &vpg) {

  std::vector<gc::Vector3> vecPosition(
      vpg.inputVertexPositions.raw().data(),
      vpg.inputVertexPositions.raw().data() +
          vpg.inputVertexPositions.raw().size());

  std::unique_ptr<gcs::SurfaceMesh> generalMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> generalVpg;

  std::tie(generalMesh, generalVpg) = gcs::makeGeneralHalfedgeAndGeometry(
      mesh.getFaceVertexList(), vecPosition);
  std::tie(L, M) = gcs::buildTuftedLaplacian(*generalMesh, *generalVpg, 1e-6);
}

/**
 * @brief Apply boundary condition mask based on boundary condition type
 *
 * @param mesh
 * @param mask
 * @param boundaryConditionType: dirichlet or neumann
 *
 */
DLL_PUBLIC inline void boundaryMask(gcs::SurfaceMesh &mesh,
                                    gcs::VertexData<bool> &mask,
                                    std::string boundaryConditionType) {
  // for (gcs::Vertex v : mesh.vertices()) {
  //   if (v.isBoundary()) {
  //     mask[v.getIndex()] = 0;
  //     for (gcs::Halfedge he : v.outgoingHalfedges()) {
  //       mask[he.next().vertex().getIndex()] = 0;
  //     }
  //   }
  // }
  if (boundaryConditionType == "fixed") {
    for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
      for (gcs::Vertex v0 : bl.adjacentVertices()) {
        for (gcs::Vertex v01 : v0.adjacentVertices()) {
          for (gcs::Vertex v012 : v01.adjacentVertices()) {
            mask[v012] = false;
          }
        }
      }
    }
  } else if (boundaryConditionType == "pin") {
    for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
      for (gcs::Vertex v0 : bl.adjacentVertices()) {
        mask[v0] = false;
      }
    }
  } else if (boundaryConditionType == "roller" ||
             boundaryConditionType == "none") {
  } else {
    throw std::runtime_error(
        "boundaryMask(bool): boundaryConditionType not defined!");
  }

  if (mask.raw().all() && boundaryConditionType != "roller" &&
      boundaryConditionType != "none") {
    std::cout << "\nboundaryMask(bool): WARNING: there is no boundary vertex "
                 "in the mesh!"
              << std::endl;
  }
  if (!mask.raw().any()) {
    std::cout << "\nboundaryMask(bool): WARNING: there is no non-masked DOF in "
                 "the mesh!"
              << std::endl;
  }
}
DLL_PUBLIC inline void boundaryMask(gcs::SurfaceMesh &mesh,
                                    gcs::VertexData<gc::Vector3> &mask,
                                    std::string boundaryConditionType) {
  if (boundaryConditionType == "fixed") {
    for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
      for (gcs::Vertex v0 : bl.adjacentVertices()) {
        for (gcs::Vertex v01 : v0.adjacentVertices()) {
          for (gcs::Vertex v012 : v01.adjacentVertices()) {
            mask[v012] = gc::Vector3{0.0, 0.0, 0.0};
          }
        }
      }
    }
  } else if (boundaryConditionType == "pin") {
    for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
      for (gcs::Vertex v0 : bl.adjacentVertices()) {
        mask[v0] = gc::Vector3{0.0, 0.0, 0.0};
      }
    }
  } else if (boundaryConditionType == "roller") {
    for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
      for (gcs::Vertex v0 : bl.adjacentVertices()) {
        mask[v0] = gc::Vector3{1.0, 1.0, 0};
      }
    }
  } else if (boundaryConditionType == "none") {
  } else {
    throw std::runtime_error(
        "boundaryMask(double): boundaryConditionType not defined!");
  }
  if (!(gc::EigenMap<double, 3>(mask).array() < 0.5).any() &&
      boundaryConditionType != "none") {
    std::cout << "\nboundaryMask(double): WARNING: there is no boundary vertex "
                 "in the mesh!"
              << std::endl;
  }
  if (!(gc::EigenMap<double, 3>(mask).array() > 0.5).any()) {
    std::cout << "\nboundaryMask(double): WARNING: there is no non-masked DOF "
                 "in the mesh!"
              << std::endl;
  }
}

/**
 * @brief Remove the rigid body translation
 *
 * @param Eigen force matrix
 */
DLL_PUBLIC inline void
removeTranslation(Eigen::Matrix<double, Eigen::Dynamic, 3> &force) {
  force = force.rowwise() - ((force).colwise().sum() / force.rows());
}
DLL_PUBLIC inline void
removeTranslation(Eigen::Matrix<double, Eigen::Dynamic, 3> &&force) {
  force = force.rowwise() - ((force).colwise().sum() / force.rows());
}

/**
 * @brief Remove the rigid body Rotation
 *
 * @param Eigen force matrix
 * @param Eigen position matrix
 */
DLL_PUBLIC inline void
removeRotation(Eigen::Matrix<double, Eigen::Dynamic, 3> position,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &force) {
  force = force.rowwise() -
          (rowwiseCrossProduct(position, force).colwise().sum() / force.rows());
}
DLL_PUBLIC inline void
removeRotation(Eigen::Matrix<double, Eigen::Dynamic, 3> position,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &&force) {
  force = force.rowwise() -
          (rowwiseCrossProduct(position, force).colwise().sum() / force.rows());
}

/**
 * @brief Gaussian distribution
 *
 * @param distance vector
 * @param standard deviation
 */
DLL_PUBLIC inline void
gaussianDistribution(Eigen::Matrix<double, Eigen::Dynamic, 1> &distribution,
                     Eigen::Matrix<double, Eigen::Dynamic, 1> distance,
                     double stdDev) {
  distribution =
      (-distance.array() * distance.array() / (2 * stdDev * stdDev)).exp() /
      (stdDev * pow(constants::PI * 2, 0.5));
}

/**
 * @brief find the closest point index to a given point
 *
 * @param mesh mesh
 * @param vpg geometry
 * @param position position of the target space point
 * @param geodesicDistance geodesic distance from a particular point in order to
 * specify range of search
 * @param range range of search
 */
DLL_PUBLIC inline gcs::Vertex
closestVertexToPt(gcs::SurfaceMesh &mesh, gcs::VertexPositionGeometry &vpg,
                  std::vector<double> position,
                  gcs::VertexData<double> &geodesicDistance,
                  double range = 1e10) {
  gcs::Vertex theVertex;
  double shorestDistance = 1e18;
  bool isIntialized = !geodesicDistance.raw().isZero(0);
  if (!isIntialized) {
    std::cout << "\nWARNING: closestVertexToPt: geodesicDistance not "
                 "initialized, searching for all "
                 "vertices!"
              << std::endl;
  }
  for (gcs::Vertex v : mesh.vertices()) {
    if (geodesicDistance[v] > range) {
      continue;
    }
    if (geodesicDistance[v] <= 0 && isIntialized) {
      std::cout
          << "\nWARNING: closestVertexToPt: geodesicDistance <= 0, may be "
             "uninitialized/updated!"
          << std::endl;
    }
    double distance;
    if (position.size() == 2) {
      distance = (gc::Vector2{vpg.inputVertexPositions[v].x,
                              vpg.inputVertexPositions[v].y} -
                  gc::Vector2{position[0], position[1]})
                     .norm();
    } else if (position.size() == 3) {
      distance = (vpg.inputVertexPositions[v] -
                  gc::Vector3{position[0], position[1], position[2]})
                     .norm();
    } else {
      throw std::runtime_error(
          "closestVertexToPt: does not support non-2d/3d position vector!");
    }
    if (distance < shorestDistance) {
      shorestDistance = distance;
      theVertex = v;
    }
  }
  return theVertex;
}

/**
 * @brief Calculate the line tension
 *
 * @param vertexPositionGeometry
 * @param mean curvature H
 * @param vertex v
 * @param lineTensionPressure
 */
DLL_PUBLIC inline void
findVertexLineTension(gcs::VertexPositionGeometry &vpg, double eta,
                      Eigen::Matrix<double, Eigen::Dynamic, 1> &H,
                      gcs::Vertex v, gcs::Halfedge isoHe, gc::Vector3 gradVec,
                      gcs::VertexData<gc::Vector3> &lineTensionPressure) {}

/**
 * @brief height = 1 tanh step function with radius r
 *
 * @param (double) sharpness of transition
 * @param (double) radius of height = 1
 * @param (Eigen vector) distance vector
 * @param (vertexPositionGeometry) vpg
 *
 */
DLL_PUBLIC inline void
tanhDistribution(gcs::VertexPositionGeometry &vpg,
                 Eigen::Matrix<double, Eigen::Dynamic, 1> &distribution,
                 Eigen::Matrix<double, Eigen::Dynamic, 1> &distance,
                 double sharpness, std::vector<double> &axes) {
  distribution.resize(distance.rows(), 1);
  if (axes[0] == axes[1]) {
    Eigen::MatrixXd radius_vec =
        Eigen::MatrixXd::Constant(distance.rows(), 1, axes[0]);
    distribution =
        0.5 *
        (1.0 + (sharpness * (radius_vec - distance)).array().tanh()).matrix();
  } else {
    double x, y, cos_t, radius;
    for (size_t i = 0; i < distance.rows(); i++) {
      x = vpg.inputVertexPositions[i].x;
      y = vpg.inputVertexPositions[i].y;
      cos_t = vpg.inputVertexPositions[i].x / sqrt(x * x + y * y);
      radius = axes[0] * axes[1] /
               sqrt((axes[0] * axes[0] - axes[1] * axes[1]) * cos_t * cos_t +
                    axes[1] * axes[1]);
      distribution[i] = 0.5 * (1 + tanh(sharpness * (radius - distance[i])));
    }
  }
}

/**
 * @brief height = 1 for elliptical domain
 *
 * @param (double) sharpness of transition
 * @param (double) radius of height = 1
 * @param (Eigen vector) distance vector
 * @param (vertexPositionGeometry) vpg
 *
 */
DLL_PUBLIC inline void
ellipticDistribution(gcs::VertexPositionGeometry &vpg,
                     Eigen::Matrix<double, Eigen::Dynamic, 1> &distribution,
                     Eigen::Matrix<double, Eigen::Dynamic, 1> &distance,
                     std::vector<double> &axes) {
  distribution.resize(distance.rows(), 1);
  if (axes[0] == axes[1]) {
    distribution = (distance.array() < axes[0]).cast<double>();
  } else {
    double x, y, cos_t, radius;
    for (size_t i = 0; i < distance.rows(); i++) {
      x = vpg.inputVertexPositions[i].x;
      y = vpg.inputVertexPositions[i].y;
      cos_t = vpg.inputVertexPositions[i].x / sqrt(x * x + y * y);
      radius = axes[0] * axes[1] /
               sqrt((axes[0] * axes[0] - axes[1] * axes[1]) * cos_t * cos_t +
                    axes[1] * axes[1]);
      distribution[i] = (double)(distance[i] < radius);
    }
  }
}

#ifdef __linux__
/**
 * @brief get the duration of world time from the start to now
 *
 * @param start reference to the starting time value
 * @return time duration
 */
DLL_PUBLIC inline double getDuration(timeval &start) {
  struct timeval end;
  gettimeofday(&end, NULL);
  return ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec -
          start.tv_usec) /
         1.e6;
}
#endif

/**
 * @brief Get the Face Vertex Matrix object
 *
 * @param mesh   Mesh of interest
 * @return Eigen matrix of uint32_t indices
 */
DLL_PUBLIC inline Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3,
                                Eigen::RowMajor>
getFaceVertexMatrix(gcs::SurfaceMesh &mesh) {
  if (!mesh.isTriangular())
    throw std::runtime_error("Mesh is not triangular.");

  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> result(
      mesh.nFaces(), 3);

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

} // namespace mem3dg
