#include "csc_matrix_conv.h"
#include <gtest/gtest.h>
#include <iostream>
#include <eigen3/Eigen/SparseCore>

class CSCMatrixConvTest : public ::testing::Test
{
public:
};

TEST_F(CSCMatrixConvTest, DenseToSparseTest1)
{
  Eigen::MatrixXd dense(2, 2);
  dense << 4, 1, 1, 2;

  osqp::CSC_Matrix sparse = osqp::calCSCMatrix(dense);

  std::vector<c_float> vals_true{4.0, 1.0, 1.0, 2.0};
  std::vector<c_int> row_idxs_true{0, 1, 0, 1};
  std::vector<c_int> col_idxs_true{0, 2, 4};

  ASSERT_TRUE(std::equal(sparse.vals.begin(), sparse.vals.end(), vals_true.begin()));
  ASSERT_TRUE(std::equal(sparse.row_idxs.begin(), sparse.row_idxs.end(), row_idxs_true.begin()));
  ASSERT_TRUE(std::equal(sparse.col_idxs.begin(), sparse.col_idxs.end(), col_idxs_true.begin()));
}

TEST_F(CSCMatrixConvTest, DenseToSparseTest2)
{
  Eigen::MatrixXd dense(4, 6);
  dense << 11, 0, 0, 14, 0, 16, 0, 22, 0, 0, 25, 26, 0, 0, 33, 34, 0, 36, 41, 0, 43, 44, 0, 46;

  osqp::CSC_Matrix sparse = osqp::calCSCMatrix(dense);

  osqp::printCSCMatrix(sparse);

  std::vector<c_float> vals_true{11, 41, 22, 33, 43, 14, 34, 44, 25, 16, 26, 36, 46};
  std::vector<c_int> row_idxs_true{0, 3, 1, 2, 3, 0, 2, 3, 1, 0, 1, 2, 3};
  std::vector<c_int> col_idxs_true{0, 2, 3, 5, 8, 9, 13};

  ASSERT_TRUE(std::equal(sparse.vals.begin(), sparse.vals.end(), vals_true.begin()));
  ASSERT_TRUE(std::equal(sparse.row_idxs.begin(), sparse.row_idxs.end(), row_idxs_true.begin()));
  ASSERT_TRUE(std::equal(sparse.col_idxs.begin(), sparse.col_idxs.end(), col_idxs_true.begin()));
}

TEST_F(CSCMatrixConvTest, TrapesoidalToSparseTest1)
{
  Eigen::MatrixXd dense(2, 2);
  dense << 4, 1, 1, 2;

  osqp::CSC_Matrix sparse = osqp::calCSCMatrixTrapesoidal(dense);

  std::vector<c_float> vals_true{4.0, 1.0, 2.0};
  std::vector<c_int> row_idxs_true{0, 0, 1};
  std::vector<c_int> col_idxs_true{0, 1, 3};

  ASSERT_TRUE(std::equal(sparse.vals.begin(), sparse.vals.end(), vals_true.begin()));
  ASSERT_TRUE(std::equal(sparse.row_idxs.begin(), sparse.row_idxs.end(), row_idxs_true.begin()));
  ASSERT_TRUE(std::equal(sparse.col_idxs.begin(), sparse.col_idxs.end(), col_idxs_true.begin()));
}

TEST_F(CSCMatrixConvTest, TrapesoidalToSparseTest2)
{
  Eigen::MatrixXd dense(4, 4);
  dense << 11, 0, 0, 14, 0, 22, 0, 0, 0, 0, 33, 34, 41, 0, 43, 44;

  osqp::CSC_Matrix sparse = osqp::calCSCMatrixTrapesoidal(dense);

  std::vector<c_float> vals_true{11, 22, 33, 14, 34, 44};
  std::vector<c_int> row_idxs_true{0, 1, 2, 0, 2, 3};
  std::vector<c_int> col_idxs_true{0, 1, 2, 3, 6};

  ASSERT_TRUE(std::equal(sparse.vals.begin(), sparse.vals.end(), vals_true.begin()));
  ASSERT_TRUE(std::equal(sparse.row_idxs.begin(), sparse.row_idxs.end(), row_idxs_true.begin()));
  ASSERT_TRUE(std::equal(sparse.col_idxs.begin(), sparse.col_idxs.end(), col_idxs_true.begin()));
}
