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

#include <array>
#include <vector>

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/halfedge_factories.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/tufted_laplacian.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"
#include "geometrycentral/utilities/vector3.h"
#include <geometrycentral/surface/surface_mesh.h>

#include <Eigen/Core>

#include <math.h>
#ifdef __linux__
#include <sys/time.h>
#endif

#include "mem3dg/constants.h"
#include "mem3dg/macros.h"
#include "mem3dg/type_utilities.h"

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
 * @brief find the range of data based on percentile (quartile)
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
 * @brief test whether exist outliers to the set of data based range function
 *
 * For example: (Outlier <--> r) < [threshold * (l <--> r)]
 *
 * @param vec data vector
 * @param threshold coefficient used to bound the outlier.
 * @return mask
 */
DLL_PUBLIC inline EigenVectorX1_T<bool> outlierMask(const EigenVectorX1d &vec,
                                                    double threshold = 0.5,
                                                    bool negate = false) {
  EigenVectorX1_T<bool> mask(vec.rows(), 1);
  EigenVectorX1d sorted_vec(vec.rows(), 1);
  sortVector(vec, sorted_vec);
  double r, l, range;
  findRange(sorted_vec.data(), sorted_vec.size(), r, l);
  range = r - l;
  if (negate) {
    // outlier: true
    for (Eigen::Index i = 0; i < vec.rows(); i++) {
      mask[i] = ((vec[i] - r < threshold * range) &&
                 (l - vec[i] < threshold * range));
    }
  } else {
    // outlier: false
    for (Eigen::Index i = 0; i < vec.rows(); i++) {
      mask[i] = ((vec[i] - r > threshold * range) ||
                 (l - vec[i] > threshold * range));
    }
  }
  return mask;
}

/**
 * @brief Signal handler for pybind
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
 * @return double
 */
DLL_PUBLIC inline double signedVolumeFromFace(gc::Vector3 &p0, gc::Vector3 &p1,
                                              gc::Vector3 &p2) {
  double v321 = p2.x * p1.y * p0.z;
  double v231 = p1.x * p2.y * p0.z;
  double v312 = p2.x * p0.y * p1.z;
  double v132 = p0.x * p2.y * p1.z;
  double v213 = p1.x * p0.y * p2.z;
  double v123 = p0.x * p1.y * p2.z;
  return (-v321 + v231 + v312 - v132 - v213 + v123) / 6.0;
}

DLL_PUBLIC inline double
signedVolumeFromFace(gcs::Vertex &&v0, gcs::Vertex &&v1, gcs::Vertex &v2,
                     gcs::VertexPositionGeometry &vpg) {
  return signedVolumeFromFace(vpg.inputVertexPositions[v0],
                              vpg.inputVertexPositions[v1],
                              vpg.inputVertexPositions[v2]);
}

DLL_PUBLIC
inline double signedVolumeFromFace(gcs::Face &f,
                                   gcs::VertexPositionGeometry &vpg) {
  gc::Vector3 p[3];
  std::size_t i = 0;
  for (gcs::Vertex v : f.adjacentVertices()) {
    p[i] = vpg.inputVertexPositions[v];
    i++;
  }
  return signedVolumeFromFace(p[0], p[1], p[2]);
}

/**
 * @brief Get the volume of a mesh object
 *
 * @param mesh Manifold surface mesh topology
 * @param vpg The vertex position geometry of the mesh
 * @param isFillHole Flag for filling holes in open mesh
 * @return The signed volume of the geometry
 */
DLL_PUBLIC inline double getMeshVolume(gcs::ManifoldSurfaceMesh &mesh,
                                       gcs::VertexPositionGeometry &vpg,
                                       bool isFillHole = true) {
  double volume = 0;

  // Throw error or fill hole for open mesh
  if (mesh.hasBoundary()) {
    if (isFillHole) {
      for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
        gcs::Vertex theVertex = bl.halfedge().tailVertex();
        for (gcs::Halfedge e : bl.adjacentHalfedges()) {
          if (e.tailVertex() != theVertex && e.tipVertex() != theVertex) {
            volume += signedVolumeFromFace(e.tailVertex(), e.tipVertex(),
                                           theVertex, vpg);
          }
        }
      }
    } else {
      mem3dg_runtime_error("Mesh is opened, not able to ",
                           "compute enclosed volume unless filled holes!");
    }
  }

  for (gcs::Face f : mesh.faces()) {
    volume += signedVolumeFromFace(f, vpg);
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
 * @brief Convert cartesian coordinates to barycentric coordinates
 *
 * @param v1 The first basis
 * @param v2 The second basis
 * @param v3 The third basis
 * @param v
 * @return DLL_PUBLIC
 */
DLL_PUBLIC inline gc::Vector3 cartesianToBarycentric(gc::Vector2 &v1,
                                                     gc::Vector2 &v2,
                                                     gc::Vector2 &v3,
                                                     gc::Vector2 &v) {
  gc::Vector3 lambda{
      ((v2.y - v3.y) * (v.x - v3.x) + (v3.x - v2.x) * (v.y - v3.y)) /
          ((v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y)), // X
      ((v3.y - v1.y) * (v.x - v3.x) + (v1.x - v3.x) * (v.y - v3.y)) /
          ((v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y)), // Y
      0};
  lambda.z = 1 - lambda.x - lambda.y;

  return lambda;
}

/**
 * @brief Get corresponding barycentric coordinate to a face
 *
 * @param baryCoords reference to Barycentric coordinate
 * @param firstHalfedge reference to the halfedge associated with the first
 * vertex of the Barycentric coordinate
 * @return
 */
DLL_PUBLIC inline gc::Vector3
correspondBarycentricCoordinates(const gc::Vector3 baryCoords_,
                                 const gcs::Halfedge &firstHalfedge) {
  std::size_t vertexInd = 0;
  gc::Vector3 baryCoords;
  for (gcs::Vertex v : firstHalfedge.face().adjacentVertices()) {
    if (v == firstHalfedge.vertex()) {
      baryCoords[vertexInd] = baryCoords_.x;
    } else if (v == firstHalfedge.next().vertex()) {
      baryCoords[vertexInd] = baryCoords_.y;
    } else if (v == firstHalfedge.next().next().vertex()) {
      baryCoords[vertexInd] = baryCoords_.z;
    } else if (firstHalfedge.face().isBoundaryLoop()) {
      mem3dg_runtime_error("Face is on boundary!");
    } else {
      mem3dg_runtime_error("Undefined behavior");
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

  // could be used to project onto other direction if needed
  return 0.5 * signedArea.norm();
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
DLL_PUBLIC inline void boundaryProteinMask(gcs::SurfaceMesh &mesh,
                                           gcs::VertexData<double> &mask,
                                           std::string boundaryConditionType) {
  // for (gcs::Vertex v : mesh.vertices()) {
  //   if (v.isBoundary()) {
  //     mask[v.getIndex()] = 0;
  //     for (gcs::Halfedge he : v.outgoingHalfedges()) {
  //       mask[he.next().vertex().getIndex()] = 0;
  //     }
  //   }
  // }
  if (boundaryConditionType == "pin") {
    for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
      for (gcs::Vertex v0 : bl.adjacentVertices()) {
        mask[v0] = 0;
      }
    }
  } else if (boundaryConditionType == "none") {
  } else {
    mem3dg_runtime_error("boundaryConditionType not defined!");
  }
  if (!(mask.raw().array() < 0.5).any() && boundaryConditionType != "none") {
    std::cout << "\nboundaryProteinMask: WARNING: there is no boundary vertex "
                 "in the mesh!"
              << std::endl;
  }
  if (!(mask.raw().array() > 0.5).any()) {
    std::cout << "\nboundaryProteinMask: WARNING: there is no non-masked DOF "
                 "in the mesh!"
              << std::endl;
  }
}
DLL_PUBLIC inline void boundaryForceMask(gcs::SurfaceMesh &mesh,
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
    mem3dg_runtime_error("boundaryConditionType is not defined");
  }
  if (!(gc::EigenMap<double, 3>(mask).array() < 0.5).any() &&
      boundaryConditionType != "none") {
    mem3dg_runtime_warning("there is no boundary vertex in the mesh!");
  }
  if (!(gc::EigenMap<double, 3>(mask).array() > 0.5).any()) {
    mem3dg_runtime_warning("there is no non-masked DOF in the mesh!");
  }
}

/**
 * @brief Remove rigid body translation
 *
 * @param force   Matrix of forces to process
 */
DLL_PUBLIC inline void removeTranslation(Eigen::Ref<EigenVectorX3dr> force) {
  force = force.rowwise() - ((force).colwise().sum() / force.rows());
}

/**
 * @brief Remove the rigid body Rotation
 *
 * @param Eigen force matrix
 * @param Eigen position matrix
 */
DLL_PUBLIC inline void
removeRotation(const Eigen::Ref<const EigenVectorX3dr> &position,
               Eigen::Ref<EigenVectorX3dr> force) {
  Eigen::Matrix<double, 1, 3> sum(0, 0, 0);
  for (Eigen::Index i = 0; i < force.rows(); ++i) {
    sum += position.row(i).cross(force.row(i));
  }
  force = force.rowwise() - (sum / force.rows());
}

/**
 * @brief Gaussian distribution
 *
 * @param distance vector
 * @param standard deviation
 */
inline double gaussianDistribution(const double &distance,
                                   const double &stdDev) {
  return exp(-distance * distance / (2 * stdDev * stdDev));
}
inline double
gaussianDistribution(const double &distance,
                     const gc::Vector3 &vertexPositionFromPtInd,
                     const std::array<gc::Vector3, 2> &tangentBasis,
                     const std::array<double, 2> &stdDev) {
  double x = gc::dot(vertexPositionFromPtInd, tangentBasis[0]);
  double y = gc::dot(vertexPositionFromPtInd, tangentBasis[1]);
  double cos_t = 1.0;
  if (!(x == 0 && y == 0)) {
    cos_t = x / sqrt(x * x + y * y);
  }
  return exp(
      -((distance * distance * cos_t * cos_t / (2 * stdDev[0] * stdDev[0])) +
        (distance * distance * (1 - cos_t * cos_t) /
         (2 * stdDev[1] * stdDev[1]))));
}
DLL_PUBLIC inline void
gaussianDistribution(Eigen::Ref<EigenVectorX1d> distribution,
                     const Eigen::Ref<const EigenVectorX1d> &distances,
                     const double stdDev) {
  distribution =
      (-distances.array() * distances.array() / (2 * stdDev * stdDev)).exp();
}
DLL_PUBLIC inline void gaussianDistribution(
    EigenVectorX1d &distribution,
    const Eigen::Ref<const EigenVectorX1d> &distances,
    const gcs::VertexData<gc::Vector3> &vertexPositionsFromPtInd,
    const std::array<gc::Vector3, 2> &tangentBasis,
    const std::array<double, 2> &stdDev) {
  distribution.resize(distances.rows(), 1);
  distribution.setConstant(1.0);
  for (Eigen::Index i = 0; i < distances.rows(); i++) {
    if (distances[i] != 0) {
      distribution[i] = gaussianDistribution(
          distances[i], vertexPositionsFromPtInd[i], tangentBasis, stdDev);
    }
  }
}

/**
 * @brief height = 1 tanh step function with radius r
 *
 * @param sharpness sharpness of transition
 * @param ax radius of height = 1
 * @param distance distance vector
 *
 */
inline double tanhDistribution(const double &distance, const double &sharpness,
                               const double &ax) {
  return 0.5 * (1.0 + tanh(sharpness * (ax - distance)));
}
inline double tanhDistribution(const double &distance,
                               const gc::Vector3 &vertexPositionFromPtInd,
                               const std::array<gc::Vector3, 2> &tangentBasis,
                               const double &sharpness,
                               const std::array<double, 2> &axes) {
  double x = gc::dot(vertexPositionFromPtInd, tangentBasis[0]);
  double y = gc::dot(vertexPositionFromPtInd, tangentBasis[1]);
  double cos_t = 1.0;
  if (!(x == 0 && y == 0)) {
    cos_t = x / sqrt(x * x + y * y);
  }
  double ax = axes[0] * axes[1] /
              sqrt((axes[0] * axes[0] - axes[1] * axes[1]) * cos_t * cos_t +
                   axes[1] * axes[1]);
  return 0.5 * (1.0 + tanh(sharpness * (ax - distance)));
}
DLL_PUBLIC inline void
tanhDistribution(Eigen::Ref<EigenVectorX1d> distribution,
                 const Eigen::Ref<const EigenVectorX1d> &distances,
                 const double sharpness, const double &ax) {
  Eigen::MatrixXd radius_vec =
      Eigen::MatrixXd::Constant(distances.rows(), 1, ax);
  distribution =
      0.5 *
      (1.0 + (sharpness * (radius_vec - distances)).array().tanh()).matrix();
}
DLL_PUBLIC inline void
tanhDistribution(EigenVectorX1d &distribution,
                 const Eigen::Ref<const EigenVectorX1d> &distances,
                 const gcs::VertexData<gc::Vector3> &vertexPositionsFromPtInd,
                 const std::array<gc::Vector3, 2> &tangentBasis,
                 const double sharpness, const std::array<double, 2> &axes) {
  distribution.resize(distances.rows(), 1);
  distribution.setConstant(1.0);
  for (Eigen::Index i = 0; i < distances.rows(); i++) {
    if (distances[i] != 0) {
      distribution[i] =
          tanhDistribution(distances[i], vertexPositionsFromPtInd[i],
                           tangentBasis, sharpness, axes);
    }
  }
}

/**
 * @brief Slice string using two deliminator
 *
 */
DLL_PUBLIC inline std::string
sliceString(std::string fileName, std::string delim1, std::string delim2) {
  // int start = 0;
  int start = fileName.find_last_of(delim1) + 1;
  int end = fileName.find_last_of(delim2);
  std::string string;
  if (start != -1 && end != -1 && end > start) {
    string = fileName.substr(start, end - start);
  }
  return string;
}

/**
 * @brief Mark the file name
 *
 * @param dirPath path of the directory
 * @param file name of the file, for example in the form of "/traj.nc"
 * @param marker_str marker used to mark the file, such as marker = "_failed"
 * results in new file name of "/traj_failed.nc"
 */
DLL_PUBLIC inline void markFileName(std::string filePath,
                                    std::string marker_str,
                                    std::string delimiter = ".") {
  size_t pos = 0;
  std::string token;
  if ((pos = filePath.find_last_of(delimiter)) != std::string::npos) {
    token = filePath.substr(0, pos);
    filePath.erase(0, pos + delimiter.length());
  }

  std::string newFilePath = token, oldFilePath = token;
  oldFilePath.append(delimiter);
  oldFilePath.append(filePath);
  newFilePath.append(marker_str);
  newFilePath.append(delimiter);
  newFilePath.append(filePath);

  // rename file
  int result = rename(oldFilePath.c_str(), newFilePath.c_str());
  if (result != 0)
    mem3dg_runtime_error("Error renaming the file!");
}

/**
 * @brief height = 1 for jump step domain
 *
 * @param (double) sharpness of transition
 * @param (double) radius of height = 1
 * @param (Eigen vector) distance vector
 * @param (vertexPositionGeometry) vpg
 *
 */
DLL_PUBLIC inline double
jumpDistribution(const double &distance,
                 const gc::Vector3 &vertexPositionFromPtInd,
                 const std::array<gc::Vector3, 2> &tangentBasis,
                 const std::array<double, 2> &axes) {
  double x = gc::dot(vertexPositionFromPtInd, tangentBasis[0]);
  double y = gc::dot(vertexPositionFromPtInd, tangentBasis[1]);
  double cos_t = 1.0;
  if (!(x == 0 && y == 0)) {
    cos_t = x / sqrt(x * x + y * y);
  }
  double ax = axes[0] * axes[1] /
              sqrt((axes[0] * axes[0] - axes[1] * axes[1]) * cos_t * cos_t +
                   axes[1] * axes[1]);
  return (double)(distance < ax);
}
DLL_PUBLIC inline void
jumpDistribution(Eigen::Ref<EigenVectorX1d> distribution,
                 const Eigen::Ref<const EigenVectorX1d> &distances,
                 const double &ax) {
  distribution = (distances.array() < ax).cast<double>();
}
DLL_PUBLIC inline void
jumpDistribution(EigenVectorX1d &distribution,
                 const Eigen::Ref<const EigenVectorX1d> &distances,
                 const gcs::VertexData<gc::Vector3> &vertexPositionsFromPtInd,
                 const std::array<gc::Vector3, 2> &tangentBasis,
                 const std::array<double, 2> &axes) {
  distribution.resize(distances.rows(), 1);
  distribution.setConstant(1.0);
  for (Eigen::Index i = 0; i < distances.rows(); i++) {
    if (distances[i] != 0) {
      distribution[i] = jumpDistribution(
          distances[i], vertexPositionsFromPtInd[i], tangentBasis, axes);
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

} // namespace mem3dg
