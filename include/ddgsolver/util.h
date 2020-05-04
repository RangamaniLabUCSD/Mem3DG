
#pragma once

#include <Eigen/Core>

#include <geometrycentral/surface/halfedge_containers.h>

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

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
 * @brief Access the raw buffer from a MeshData object
 *
 * @tparam E    Typename of the mesh element
 * @tparam T    Typename of the stored data
 * @param data  The data
 * @return T*   C-style array to raw buffer
 */
template <typename E, typename T>
inline T *getBuffer(gcs::MeshData<E, T> &data) {
  return data.rawdata().data();
}

/**
 * @brief Access the raw buffer from a MeshData object
 *
 * @tparam E    Typename of the mesh element
 * @tparam T    Typename of the stored data
 * @param data  The data
 * @return T*   C-style array to raw buffer
 */
template <typename E, typename T>
inline const T *getBuffer(const gcs::MeshData<E, T> &data) {
  return data.rawdata().data();
}

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
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T        Typename of the output data
 * @tparam k        Number of columns in the output
 * @tparam E        Typename of the mesh element
 * @tparam O        Typename of the input data
 * @tparam Options  Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @param vec       The data
 * @return AlignedEigenMap_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename E, typename O,
          int Options = Eigen::RowMajor>
AlignedEigenMap_T<T, k, Options> EigenMap(gcs::MeshData<E, O> &vec) {
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");
  return AlignedEigenMap_T<T, k, Options>(reinterpret_cast<T *>(getBuffer(vec)),
                                          vec.size(), k);
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T        Typename of the output data
 * @tparam k        Number of columns in the output
 * @tparam E        Typename of the mesh element
 * @tparam O        Typename of the input data
 * @tparam Options  Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @param vec       The data
 * @return AlignedEigenMap_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename E, typename O,
          int Options = Eigen::RowMajor>
ConstAlignedEigenMap_T<T, k, Options> EigenMap(const gcs::MeshData<E, O> &vec) {
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");
  return ConstAlignedEigenMap_T<T, k, Options>(
      reinterpret_cast<const T *>(getBuffer(vec)), vec.size(), k);
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
 * @tparam T    Typename of the stored data
 * @param vec   The data
 * @return AlignedEigenMap_T<T, 1>    Eigen map object to the buffer
 */
template <typename E, typename T>
AlignedEigenMap_T<T, 1> EigenMap(gcs::MeshData<E, T> &vec) {
  return AlignedEigenMap_T<T, 1>(getBuffer(vec), vec.size());
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * @tparam T    Typename of the stored data
 * @param vec   The data
 * @return AlignedEigenMap_T<T, 1>    Eigen map object to the buffer
 */
template <typename E, typename T>
ConstAlignedEigenMap_T<T, 1> EigenMap(const gcs::MeshData<E, T> &vec) {
  return ConstAlignedEigenMap_T<T, 1>(getBuffer(vec), vec.size());
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
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T    Typename of the output data
 * @tparam k    Number of columns in the output
 * @tparam E    Typename of the mesh element
 * @tparam O    Typename of the input data
 * @param vec   The data
 * @return AlignedEigenMap_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename E, typename O>
AlignedEigenMap_T<T, 1> FlattenedEigenMap(gcs::MeshData<E, O> &vec) {
  // Imperfect check for padding
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");
  return AlignedEigenMap_T<T, 1>(reinterpret_cast<T *>(getBuffer(vec)),
                                 k * vec.size());
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T    Typename of the output data
 * @tparam k    Number of columns in the output
 * @tparam E    Typename of the mesh element
 * @tparam O    Typename of the input data
 * @param vec   The data
 * @return AlignedEigenMap_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename E, typename O>
ConstAlignedEigenMap_T<T, 1> FlattenedEigenMap(const gcs::MeshData<E, O> &vec) {
  // Imperfect check for padding
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");
  return ConstAlignedEigenMap_T<T, 1>(reinterpret_cast<const T *>(getBuffer(vec)),
                                      k * vec.size());
}
} // end namespace ddgsolver
