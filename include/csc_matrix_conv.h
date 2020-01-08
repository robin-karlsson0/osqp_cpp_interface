#ifndef CSC_MATRIX_CONV_H
#define CSC_MATRIX_CONV_H

#include <eigen3/Eigen/Core>
#include <iostream>
#include <vector>
#include <exception>
#include "osqp.h" // for 'c_int' type ('long' or 'long long')

namespace osqp
{

// Struct for containing a 'Compressed-Column-Sparse Matrix' object.
//
// Elements:
//   vals: Vector of non-zero values. Ex: [4,1,1,2]
//   row_idxs:  Vector of row index corresponding to values. Ex: [0, 1, 0, 1] (Eigen: 'inner')
//   col_idxs:  Vector of 'val' indices where each column starts. Ex: [0, 2, 4] (Eigen: 'outer')
struct CSC_Matrix
{
  std::vector<c_float> vals;
  std::vector<c_int> row_idxs;
  std::vector<c_int> col_idxs;
};

CSC_Matrix calCSCMatrix(const Eigen::MatrixXd &mat);
CSC_Matrix calCSCMatrixTrapesoidal(const Eigen::MatrixXd &mat);

void printCSCMatrix(CSC_Matrix &csc_mat);

} // namespace osqp

#endif // CSC_MATRIX_CONV_H
