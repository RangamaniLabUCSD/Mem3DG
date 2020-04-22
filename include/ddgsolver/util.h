
#pragma once

#include <Eigen/Core>

template <typename T>
using AlignedVector_t = typename std::vector<T, Eigen::aligned_allocator<T> >;

template <typename T, std::size_t k>
using EigenMapFromAlignedVector = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, k>, Eigen::AlignedMax>;

template <typename T, std::size_t k, typename O>
EigenMapFromAlignedVector<T, k> mapVecToEigen(AlignedVector_t<O>& vec){
    // Imperfect check for padding
    static_assert(std::is_standard_layout<T>::value && std::is_trivially_copyable<O>::value, "O must be a POD type.");
    static_assert(sizeof(O) == k*sizeof(T), "sizeof(O) must be a k multiple of sizeof(T)");

    EigenMapFromAlignedVector<T, k> map(reinterpret_cast<T*>(vec.data()
        ), vec.size(), k);

    return map;
}