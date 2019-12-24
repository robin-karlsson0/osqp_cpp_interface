#include "osqp_interface.h"
#include <gtest/gtest.h>
#include <iostream>
#include <eigen3/Eigen/SparseCore>

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

class OSQPInterfaceTest : public ::testing::Test
{
public:
};

TEST_F(OSQPInterfaceTest, initializeInterfaceWithoutData)
{
  osqp::OSQPInterface osqp_interface = osqp::OSQPInterface(1e-6);
  c_int exitflag = osqp_interface.getExitFlag();

  // Expect interface to finish initialization and set flag to '0'
  ASSERT_EQ(exitflag, 0);
}

TEST_F(OSQPInterfaceTest, initiazeInterfaceWithData)
{
  // Test problem
  Eigen::MatrixXd P(2, 2);
  P << 1, -1, -1, 2;
  std::vector<double> q{-2, -6};
  std::vector<double> l{-1000, -1000, -1000};
  Eigen::MatrixXd A(3, 2);
  A << 1, 1, -1, 2, 2, 1;
  std::vector<double> u{2, 2, 3};

  osqp::OSQPInterface osqp_interface = osqp::OSQPInterface(P, A, q, l, u, 1e-6);
  c_int exitflag = osqp_interface.getExitFlag();

  // Expect interface to finish initialization and set flag to '0'
  ASSERT_EQ(exitflag, 0);
}

TEST_F(OSQPInterfaceTest, optimizeInitializedProblem)
{
  // Test problem
  Eigen::MatrixXd P(2, 2);
  P << 1, -1, -1, 2;
  std::vector<double> q{-2, -6};
  std::vector<double> l{-1000, -1000, -1000};
  Eigen::MatrixXd A(3, 2);
  A << 1, 1, -1, 2, 2, 1;
  std::vector<double> u{2, 2, 3};

  // Initialize interface
  osqp::OSQPInterface osqp_interface = osqp::OSQPInterface(P, A, q, l, u, 1e-6);

  // Optimize problem
  std::tuple<std::vector<double>, std::vector<double>> result = osqp_interface.optimize();

  // Output solution
  std::vector<double> param = std::get<0>(result);
  double x_0 = param[0];
  double x_1 = param[1];

  EXPECT_NEAR(x_0, 0.6667, 1e-4);
  EXPECT_NEAR(x_1, 1.3333, 1e-4);
}

TEST_F(OSQPInterfaceTest, optimizeUninitializedProblem)
{
  // Test problem
  Eigen::MatrixXd P(2, 2);
  P << 1, -1, -1, 2;
  std::vector<double> q{-2, -6};
  std::vector<double> l{-1000, -1000, -1000};
  Eigen::MatrixXd A(3, 2);
  A << 1, 1, -1, 2, 2, 1;
  std::vector<double> u{2, 2, 3};

  // Initialize interface
  osqp::OSQPInterface osqp_interface = osqp::OSQPInterface(1e-6);

  // Optimize problem
  std::tuple<std::vector<double>, std::vector<double>> result = osqp_interface.optimize(P, A, q, l, u);

  // Output solution
  std::vector<double> param = std::get<0>(result);
  double x_0 = param[0];
  double x_1 = param[1];

  EXPECT_NEAR(x_0, 0.6667, 1e-4);
  EXPECT_NEAR(x_1, 1.3333, 1e-4);
}

TEST_F(OSQPInterfaceTest, optimizeAfterUpdates)
{
  // Test problem
  Eigen::MatrixXd P(2, 2);
  P << 1, -1, -1, 2;
  std::vector<double> q{-2, -6};
  std::vector<double> l{-1000, -1000, -1000};
  Eigen::MatrixXd A(3, 2);
  A << 1, 1, -1, 2, 2, 1;
  std::vector<double> u{2, 2, 3};

  // Initialize interface
  osqp::OSQPInterface osqp_interface = osqp::OSQPInterface(P, A, q, l, u, 1e-6);

  // Optimize problem
  std::tuple<std::vector<double>, std::vector<double>> result = osqp_interface.optimize();

  // Update P matrix
  osqp_interface.updateP(P);
  osqp_interface.updateA(A);
  osqp_interface.updateQ(q);
  osqp_interface.updateL(l);
  osqp_interface.updateU(u);

  // Optimize problem
  result = osqp_interface.optimize();

  // Output solution
  std::vector<double> param = std::get<0>(result);
  double x_0 = param[0];
  double x_1 = param[1];

  EXPECT_NEAR(x_0, 0.6667, 1e-4);
  EXPECT_NEAR(x_1, 1.3333, 1e-4);
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
