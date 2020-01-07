#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Cholesky>
#include <vector>
#include <tuple>

#include "osqp_interface.h"

// TEST PROBLEM
//  ref: https://www.mathworks.com/help/optim/ug/quadprog.html
//
// 'P' matrix:
//   | 1  -1 |
//   | -1  2 |
//
// 'q' vector:
//   | -2 |
//   | -6 |
//
// 'l' vector:
//   | -1000 |
//   | -1000 |
//
// 'A' vector:
//   |  1  1 |
//   | -1  2 |
//   |  2  1 |
//
// 'u' vector:
//   | 2 |
//   | 2 |
//   | 3 |
//
// Solution: [0.6666..., 1.3333...]

int main(int argc, char *argv[])
{
  // Test problem
  Eigen::MatrixXd P(2, 2);
  P << 1, -1, -1, 2;
  std::vector<double> q{-2, -6};
  std::vector<double> l{-1000, -1000, -1000};
  Eigen::MatrixXd A(3, 2);
  A << 1, 1, -1, 2, 2, 1;
  std::vector<double> u{2, 2, 3};

  // Initialize the interface
  float eps_abs = 1e-6;
  osqp::OSQPInterface osqp_interface = osqp::OSQPInterface(P, A, q, l, u, eps_abs);

  // Run the optimizer
  std::tuple<std::vector<double>, std::vector<double>, int> result = osqp_interface.optimize();

  // Output solution
  std::vector<double> param = std::get<0>(result);
  for (std::vector<double>::const_iterator it = param.begin(); it != param.end(); it++)
  {
    std::cout << *it << std::endl;
  }
  int status_polish = std::get<2>(result);
  std::cout << "status_polish: " << status_polish << " (1: successful, 0: unperformed, -1: unsuccessful)" << std::endl;

  return 0;
}