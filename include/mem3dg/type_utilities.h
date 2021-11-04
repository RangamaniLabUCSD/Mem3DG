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

/**
 * @file util.h
 * @brief Utility functions primarily for interfacing with Eigen and raw buffers
 *
 */

#pragma once

#include <Eigen/Core>

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/macros.h"

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

using EigenVectorX1d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
template <typename T>
using EigenVectorX1_T = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using EigenVectorX1i = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using EigenVectorX3dr =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3ur =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3u = Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3>;

template <typename T, int k>
using EigenVectorXkr_T = Eigen::Matrix<T, Eigen::Dynamic, k, Eigen::RowMajor>;

/// Type alias for aligned vectors
template <typename T>
using AlignedVector_T = std::vector<T, Eigen::aligned_allocator<T>>;

/**
 * @brief Typename for an Eigen map to an aligned raw buffer
 *
 * This assumes that the underlying buffer is aligned with Eigen::AlignedMax.
 * Unfortunately it doesn't seem that the compiler can infer the alignment of
 * the data for you.
 *
 * @tparam T          Typename of the contained data
 * @tparam k          Number of columns
 * @tparam Options    Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 */
template <typename T, std::size_t k, int Options = Eigen::ColMajor>
using AlignedEigenMap_T =
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, k, Options>, Eigen::AlignedMax>;

/**
 * @brief Typename for an Eigen map to a const aligned raw buffer
 *
 * This assumes that the underlying buffer is aligned with Eigen::AlignedMax.
 * Unfortunately it doesn't seem that the compiler can infer the alignment of
 * the data for you.
 *
 * @tparam T        Typename of the contained data
 * @tparam k        Numver of columns
 * @tparam Options  Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 */
template <typename T, std::size_t k, int Options = Eigen::ColMajor>
using ConstAlignedEigenMap_T =
    Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, k, Options>,
               Eigen::AlignedMax>;

/**
 * @brief Compute the rowwise dot product between two Eigen Matrices
 *
 * @tparam DerivedA  Template type value of matrix A
 * @tparam DerivedB  Template type value of matrix B
 * @param A         the first matrix
 * @param B         The other matrix
 * @return auto     Intermediate return value for Eigen optimization
 */
template <typename DerivedA, typename DerivedB>
auto rowwiseDotProduct(const Eigen::DenseBase<DerivedA> &A,
                       const Eigen::DenseBase<DerivedB> &B) {
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(DerivedA, DerivedB);
  return (A.derived().array() * B.derived().array()).rowwise().sum().matrix();
}

/**
 * @brief Compute the rowwise scalar product between a matrix and scalar
 *
 * @tparam VectorType Template type of vector type
 * @tparam Derived    Template type of matrix B
 * @param A           The vector of scalar values
 * @param B           The vector of matrix values
 * @return auto
 */
template <typename VectorType, typename Derived>
auto rowwiseScalarProduct(const Eigen::DenseBase<VectorType> &A,
                          const Eigen::DenseBase<Derived> &B) {
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorType);
  if (A.rows() != B.rows()) {
    mem3dg_runtime_error("Mismatched row size, ", A.rows(), " rows in vector ",
                         B.rows(), " in matrix, for rowwise scalar product");
  }
  return (B.derived().array().colwise() * A.derived().array()).matrix();
}

/**
 * @brief Compute the rowwise cross product between two matrices
 *
 * Note that this function does not return a result of an intermediate Eigen
 * operation. Owing to some limitations, we return an evaluated temporary.
 *
 *
 * @tparam Derived  Template type value of the matrices
 * @param A         Matrix A
 * @param B         Matrix B
 * @return          Matrix of rowwise cross products
 */
template <typename Derived>
typename Derived::PlainMatrix
rowwiseCrossProduct(const Eigen::DenseBase<Derived> &A,
                    const Eigen::DenseBase<Derived> &B) {
  if (A.rows() != B.rows()) {
    mem3dg_runtime_error("Mismatched rows, ", A.rows(), " rows in A", B.rows(),
                         " in B, for rowwise cross product");
  }
  typename Derived::PlainMatrix C;
  C.resize(A.rows(), 3);
  for (std::size_t i = 0; i < A.rows(); i++) {
    C.row(i) = (A.derived().matrix().row(i)).cross(B.derived().matrix().row(i));
  }
  return C;
}

// inline EigenVectorX3dr
// rowwiseCrossProduct(const Eigen::Ref<const EigenVectorX3dr> &A,
//                     const Eigen::Ref<const EigenVectorX3dr> &B) {
//   if (A.rows() != B.rows()) {
//     mem3dg_runtime_error("Mismatched rows, ", A.rows(), " rows in A",
//     B.rows(),
//                          " in B, for rowwise cross product");
//   }
//   EigenVectorX3dr C;
//   C.resize(A.rows(), 3);
//   for (std::size_t i = 0; i < A.rows(); i++) {
//     C.row(i) =
//     (A.derived().matrix().row(i)).cross(B.derived().matrix().row(i));
//   }
//   return C;
// }

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 * Reads the data in Row Major order by default since expansion of O to k*T
 * occurs by row in theory.
 *
 * @tparam T        Typename of the output data
 * @tparam k        Number of columns in the output
 * @tparam O        Typename of the input data
 * @tparam Options  Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @param vec       The data
 * @return AlignedEigenMap_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename O, int Options = Eigen::RowMajor>
AlignedEigenMap_T<T, k, Options> EigenMap(AlignedVector_T<O> &vec) {
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");
  return AlignedEigenMap_T<T, k, Options>(reinterpret_cast<T *>(vec.data()),
                                          vec.size(), k);
}

/**
 * @brief Generate an Eigen Map to a const aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 * Reads the data in Row Major order by default since expansion of O to k*T
 * occurs by row in theory.
 *
 * @tparam T        Typename of the output data
 * @tparam k        Number of columns in the output
 * @tparam O        Typename of the input data
 * @tparam Options  Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @param vec       The data
 * @return ConstAlignedEigenMap_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename O, int Options = Eigen::RowMajor>
ConstAlignedEigenMap_T<T, k, Options> EigenMap(const AlignedVector_T<O> &vec) {
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");
  return ConstAlignedEigenMap_T<T, k, Options>(
      reinterpret_cast<const T *>(vec.data()), vec.size(), k);
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * @tparam T    Typename of the output data
 * @param vec   The data
 * @return AlignedEigenMap_T<T, 1>    Eigen map object to the buffer
 */
template <typename T>
AlignedEigenMap_T<T, 1> EigenMap(AlignedVector_T<T> &vec) {
  return AlignedEigenMap_T<T, 1>(vec.data(), vec.size());
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * @tparam T    Typename of the output data
 * @param vec   The data
 * @return AlignedEigenMap_T<T, 1>    Eigen map object to the buffer
 */
template <typename T>
ConstAlignedEigenMap_T<T, 1> EigenMap(const AlignedVector_T<T> &vec) {
  return ConstAlignedEigenMap_T<T, 1>(vec.data(), vec.size());
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T    Typename of the output data
 * @tparam k    Number of columns in the output
 * @tparam O    Typename of the input data
 * @param vec   The data
 * @return AlignedEigenMap_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename O>
AlignedEigenMap_T<T, 1> FlattenedEigenMap(AlignedVector_T<O> &vec) {
  // Imperfect check for padding
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");
  return AlignedEigenMap_T<T, 1>(reinterpret_cast<T *>(vec.data()),
                                 k * vec.size());
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T    Typename of the output data
 * @tparam k    Number of columns in the output
 * @tparam O    Typename of the input data
 * @param vec   The data
 * @return AlignedEigenMap_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename O>
ConstAlignedEigenMap_T<T, 1> FlattenedEigenMap(const AlignedVector_T<O> &vec) {
  // Imperfect check for padding
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");
  return ConstAlignedEigenMap_T<T, 1>(reinterpret_cast<const T *>(vec.data()),
                                      k * vec.size());
}

/**
 * @brief Flatten a NxK row matrix into KNx1
 *
 * @tparam Derived Typename of the Eigen dense type
 * @param matrix   The data
 * @return AlignedEigenMap_T<typename Derived::Scalar, 1> Flattened result
 */
template <typename Derived>
AlignedEigenMap_T<typename Derived::Scalar, 1, Eigen::ColMajor>
flatten(Eigen::DenseBase<Derived> &matrix) {
  static_assert((Derived::PlainMatrix::Options & Eigen::RowMajor) == 1,
                "Input must be in row-major order");
  return AlignedEigenMap_T<typename Derived::Scalar, 1, Eigen::ColMajor>(
      matrix.derived().data(), matrix.derived().size());
}

/**
 * @brief  Flatten a NxK row matrix into KNx1
 *
 * @tparam Derived  Typename of the Eigen dense type
 * @param matrix    The data
 * @return ConstAlignedEigenMap_T<typename Derived::Scalar, 1, Eigen::ColMajor>
 * Flattened result
 */
template <typename Derived>
ConstAlignedEigenMap_T<typename Derived::Scalar, 1, Eigen::ColMajor>
flatten(const Eigen::DenseBase<const Derived> &matrix) {
  static_assert((Derived::PlainMatrix::Options & Eigen::RowMajor) == 1,
                "Input must be in row-major order");
  return ConstAlignedEigenMap_T<typename Derived::Scalar, 1, Eigen::ColMajor>(
      matrix.derived().data(), matrix.derived().size());
}

template <std::size_t k, typename Derived>
AlignedEigenMap_T<typename Derived::Scalar, k, Eigen::RowMajor>
unflatten(Eigen::DenseBase<Derived> &matrix) {
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
  assert(matrix.derived().size() % k == 0);
  return AlignedEigenMap_T<typename Derived::Scalar, k, Eigen::RowMajor>(
      matrix.derived().data(), matrix.derived().size() / k, k);
}

template <std::size_t k, typename Derived>
ConstAlignedEigenMap_T<typename Derived::Scalar, k, Eigen::RowMajor>
unflatten(const Eigen::DenseBase<Derived> &matrix) {
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
  assert(matrix.derived().size() % k == 0);
  return ConstAlignedEigenMap_T<typename Derived::Scalar, k, Eigen::RowMajor>(
      matrix.derived().data(), matrix.derived().size() / k, k);
}

inline auto toMatrix(gcs::VertexData<gc::Vector3> &vector) {
  return gc::EigenMap<double, 3>(vector);
}

inline auto toMatrix(gcs::VertexData<double> &vector) { return vector.raw(); }

} // namespace mem3dg
