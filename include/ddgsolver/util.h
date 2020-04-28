
#pragma once

#include <Eigen/Core>

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
 * @tparam T    Typename of the contained data
 * @tparam k    Numver of columns
 */
template <typename T, std::size_t k>
using EigenMapFromAlignedVector_T =
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, k>, Eigen::AlignedMax>;

/**
 * @brief Access the raw buffer from a MeshData object
 *
 * @tparam E    Typename of the mesh element
 * @tparam T    Typename of the stored data
 * @param data  The data
 * @return T*   C-style array to raw buffer
 */
template <typename E, typename T>
inline T *getBuffer(gcs::MeshData<E, T> &data)
{
  return data.rawdata().data();
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
 * @return EigenMapFromAlignedVector_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename O>
EigenMapFromAlignedVector_T<T, k> mapVecToEigen(AlignedVector_T<O> &vec)
{
  // Imperfect check for padding
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");

  EigenMapFromAlignedVector_T<T, k> map(reinterpret_cast<T *>(vec.data()),
                                        vec.size(), k);

  return map;
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
 * @return EigenMapFromAlignedVector_T<T, k>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename E, typename O>
EigenMapFromAlignedVector_T<T, k> mapVecToEigen(gcs::MeshData<E, O> &vec)
{
  // Imperfect check for padding
  static_assert(std::is_standard_layout<T>::value &&
                    std::is_trivially_copyable<O>::value,
                "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T),
                "sizeof(O) must be a k multiple of sizeof(T)");

  EigenMapFromAlignedVector_T<T, k> map(reinterpret_cast<T *>(getBuffer(vec)),
                                        vec.size(), k);

  return map;
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * @tparam T    Typename of the output data
 * @param vec   The data
 * @return EigenMapFromAlignedVector_T<T, 1>    Eigen map object to the buffer
 */
template <typename T>
EigenMapFromAlignedVector_T<T, 1> mapVecToEigen(AlignedVector_T<T> &vec)
{
  EigenMapFromAlignedVector_T<T, 1> map(vec.data(), vec.size());
  return map;
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * @tparam T    Typename of the stored data
 * @param vec   The data
 * @return EigenMapFromAlignedVector_T<T, 1>    Eigen map object to the buffer
 */
template <typename E, typename T>
EigenMapFromAlignedVector_T<T, 1> mapVecToEigen(gcs::MeshData<E, T> &vec)
{
  EigenMapFromAlignedVector_T<T, 1> map(getBuffer(vec), vec.size());
  return map;
}
