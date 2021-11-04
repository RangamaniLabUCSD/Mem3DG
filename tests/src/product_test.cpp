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

#include <iostream>

#include <gtest/gtest.h>

#include "mem3dg/macros.h"
#include "mem3dg/type_utilities.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class ProductTest : public ::testing::Test {
public:
  ProductTest() {
    for (int i = 0; i < rows; ++i) {
      A(i, 0) = 3 * i;
      Am(i, 0) = 3 * i;
      A(i, 1) = 3 * i + 1;
      Am(i, 1) = 3 * i + 1;
      A(i, 2) = 3 * i + 2;
      Am(i, 2) = 3 * i + 2;
      S(i) = i;
    }

    for (int i = 0; i < rows; ++i) {
      B(i, 0) = 10 * i;
      Bm(i, 0) = 10 * i;
      B(i, 1) = 10 * i + 1;
      Bm(i, 1) = 10 * i + 1;
      B(i, 2) = 10 * i + 2;
      Bm(i, 2) = 10 * i + 2;
    }
  }
  void SetUp() {}
  void TearDown() {}

  static std::size_t const rows = 10;
  Eigen::Array<double, rows, 3> A, B;
  Eigen::Matrix<double, rows, 3> Am, Bm;
  Eigen::Matrix<double, rows, 1> S;
};

TEST_F(ProductTest, RowwiseDotProductTest) {
  auto res = mem3dg::rowwiseDotProduct(A, B);
  ASSERT_EQ(A.rows(), res.rows());
  ASSERT_EQ(res.cols(), 1);

  auto resAM = mem3dg::rowwiseDotProduct(Am, B);
  // ASSERT_EQ(A.rows(), resAM.rows());
  ASSERT_EQ(resAM.cols(), 1);

  auto resMA = mem3dg::rowwiseDotProduct(A, Bm);
  ASSERT_EQ(A.rows(), resMA.rows());
  ASSERT_EQ(resMA.cols(), 1);

  auto resm = mem3dg::rowwiseDotProduct(Am, Bm);
  ASSERT_EQ(Am.rows(), resm.rows());
  ASSERT_EQ(resm.cols(), 1);
  ASSERT_EQ(true, res.isApprox(resm));
  ASSERT_EQ(true, res.isApprox(resAM));
  ASSERT_EQ(true, res.isApprox(resMA));

  Eigen::Matrix<double, rows, 1> manual;
  for (int i = 0; i < 10; ++i) {
    manual(i) = Am(i, 0) * Bm(i, 0) + Am(i, 1) * Bm(i, 1) + Am(i, 2) * Bm(i, 2);
  }
  ASSERT_EQ(true, manual.isApprox(resm));
}

TEST_F(ProductTest, RowwiseScalarProductTest) {
  // test vector--array
  auto resA = mem3dg::rowwiseScalarProduct(S, A);
  ASSERT_EQ(A.rows(), resA.rows());
  ASSERT_EQ(resA.cols(), 3);

  // test vecotr--matrix
  auto resAm = mem3dg::rowwiseScalarProduct(S, Am);
  ASSERT_EQ(Am.rows(), resAm.rows());
  ASSERT_EQ(resAm.cols(), 3);

  ASSERT_EQ(true, resA.isApprox(resAm));

  auto resB = mem3dg::rowwiseScalarProduct(S, B);
  ASSERT_EQ(B.rows(), resB.rows());
  ASSERT_EQ(resB.cols(), 3);

  auto resBm = mem3dg::rowwiseScalarProduct(S, Bm);
  ASSERT_EQ(Bm.rows(), resBm.rows());
  ASSERT_EQ(resBm.cols(), 3);

  ASSERT_EQ(true, resB.isApprox(resBm));

  // Check that the values are correct wrt to numerical benchmark
  Eigen::Matrix<double, rows, 3> manual;
  for (std::size_t i = 0; i < 10; ++i) {
    for (std::size_t j = 0; j < 3; ++j)
      manual(i, j) = Am(i, j) * S(i);
  }
  ASSERT_EQ(true, manual.isApprox(resAm));

  // Check that runtime throws an error on mismatched size
  Eigen::Matrix<double, Eigen::Dynamic, 1> S_wrongsize;
  S_wrongsize.resize(rows + 10, 1);
  ASSERT_THROW(mem3dg::rowwiseScalarProduct(S_wrongsize, A),
               std::runtime_error);
}

TEST_F(ProductTest, RowwiseCrossProductTest) {
  auto res = mem3dg::rowwiseCrossProduct(A, B);
  ASSERT_EQ(A.rows(), res.rows());
  ASSERT_EQ(res.cols(), 3);

  for (std::size_t i = 0; i < rows; ++i) {
    ASSERT_EQ(
        true,
        (res.row(i)).isApprox(A.row(i).matrix().cross(B.row(i).matrix())));
  }

  auto res2 = mem3dg::rowwiseCrossProduct(Am, Bm);
  ASSERT_EQ(A.rows(), res2.rows());
  ASSERT_EQ(res2.cols(), 3);
}
