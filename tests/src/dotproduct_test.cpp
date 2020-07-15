

#include <iostream>

#include <gtest/gtest.h>

#include "ddgsolver/util.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class DotProductTest : public ::testing::Test {};

TEST_F(DotProductTest, DotProductTest) {
  const static std::size_t sz = 10;
  Eigen::Array<double, sz, 3> A, B;
  Eigen::Matrix<double, sz, 3> Am, Bm;

  for (int i = 0; i < 10; ++i) {
    A(i, 0) = 3 * i;
    Am(i, 0) = 3 * i;
    A(i, 1) = 3 * i + 1;
    Am(i, 1) = 3 * i + 1;
    A(i, 2) = 3 * i + 2;
    Am(i, 2) = 3 * i + 2;
  }

  for (int i = 0; i < 10; ++i) {
    B(i, 0) = 10 * i;
    Bm(i, 0) = 10 * i;
    B(i, 1) = 10 * i + 1;
    Bm(i, 1) = 10 * i + 1;
    B(i, 2) = 10 * i + 2;
    Bm(i, 2) = 10 * i + 2;
  }

  auto res = ddgsolver::dot(A, B);
  ASSERT_EQ(A.rows(), res.rows());
  ASSERT_EQ(res.cols(), 1);

  auto resm = ddgsolver::dot(Am, Bm);
  ASSERT_EQ(Am.rows(), resm.rows());
  ASSERT_EQ(resm.cols(), 1);
  ASSERT_EQ(true, res.matrix().isApprox(resm));

  Eigen::Matrix<double, sz, 1> manual;
  for (int i = 0; i < 10; ++i) {
    manual(i) = Am(i, 0) * Bm(i, 0) + Am(i, 1) * Bm(i, 1) + Am(i, 2) * Bm(i, 2);
  }
  ASSERT_EQ(true, manual.isApprox(resm));
}
} // end namespace ddgsolver
