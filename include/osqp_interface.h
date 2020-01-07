#ifndef OSQP_INTERFACE_H
#define OSQP_INTERFACE_H

#include <eigen3/Eigen/SparseCore>
#include <vector>
#include <tuple>
#include "osqp.h"

namespace osqp
{

// Struct for containing a 'Compressed-Column-Sparse Matrix' object.
//
// Elements:
//   elem_val: Vector of non-zero values. Ex: [4,1,1,2]
//   row_idx:  Row index corresponding to values. Ex: [0, 1, 0, 1] (Eigen: 'inner')
//   col_idx:  List of 'val' indices where each column starts. Ex: [0, 2, 4] (Eigen: 'outer')
struct CSC_Matrix
{
  std::vector<double> elem_val;
  std::vector<c_int> row_idx;
  std::vector<c_int> col_idx;
};

/**
  * Implementation of a native C++ interface for the OSQP solver.
  * 
  * The interface takes in the problem formalation as Eigen matrices and vectors, converts these objects into C-style
  * CSC matrices and dynamic arrays, loads the data into the OSQP workspace dataholder, and runs the optimizer.
  * 
  * The optimization results are return as a vector tuple by the optimization function.
  *   std::tuple<std::vector<double>, std::vector<double>> result = osqp_interface.optimize();
  *   std::vector<double> param = std::get<0>(result);
  *   double x_0 = param[0];
  *   double x_1 = param[1];
  * 
  * The interface can be used in several ways:
  * 
  *   1. Initialize the interface WITHOUT data. Load the problem formulation at the optimization call.
  *        osqp_interface = OSQPInterface();
  *        osqp_interface.optimize(P, A, q, l, u);
  * 
  *   2. Initialize the interface WITH data.
  *        osqp_interface = OSQPInterface(P, A, q, l, u);
  *        osqp_interface.optimize();
  * 
  *   3. WARM START OPTIMIZATION by modifying the problem formulation between optimization runs.
  *        osqp_interface = OSQPInterface(P, A, q, l, u);
  *        osqp_interface.optimize();
  *        while()
  *        {
  *          osqp_interface.updateP(P_new);
  *          osqp_interface.updateA(A_new);
  *          osqp_interface.updateQ(q_new);
  *          osqp_interface.updateL(l_new);
  *          osqp_interface.updateU(u_new);
  *          osqp_interface.optimize();
  *        }
  * 
  * Ref: https://osqp.org/docs/solver/index.html
  */
class OSQPInterface
{

private:
  /*****************************
   * OSQP WORKSPACE STRUCTURES
   *****************************/
  OSQPWorkspace *work;
  OSQPSettings *settings;
  OSQPData *data;
  // Number of parameters to optimize
  c_int param_n;

  // For destructor to know if matrices P, A are in
  bool problem_in_memory = false;

  // Runs the solver on the stored problem.
  std::tuple<std::vector<double>, std::vector<double>, int> solve();

  /*****************************
   * DATA CONVERSION FUNCTIONS
   *****************************/
  // Converts problem input matrices to CSC matrix structs.
  CSC_Matrix transformP(const Eigen::MatrixXd &P, int *nonzeros);
  CSC_Matrix transformA(const Eigen::MatrixXd &A);
  // Converts problem input vectors to dynamic arrays.
  double *transformQ(const std::vector<double> &q);
  double *transformL(const std::vector<double> &l);
  double *transformU(const std::vector<double> &u);
  // Converts an Eigen matrix into a CSC matrix struct.
  CSC_Matrix convEigenMatrixToCSCMatrix(const Eigen::MatrixXd A);
  // Converts an Eigen vector matrix into a dynamic array.
  double *convEigenVecToDynFloatArray(const Eigen::MatrixXd x);

  // Exitflag
  c_int exitflag;

  inline bool isEqual(double x, double y);

public:
  // Returns a flag for asserting interface condition (Healthy condition: 0).
  c_int getExitFlag(void);

  /****************************
   * INITIALIZATION FUNCTIONS
   ****************************/

  // Initializes the OSQP interface without setting up the problem.
  //
  // Steps:
  //   1. Initializes the OSQP object (incl. settings, data objects).
  //   2. Solver settings (accuracy etc.).
  OSQPInterface(const c_float eps_abs);

  // Initializes the OSQP solver interface and sets up the problem.
  //
  // Steps:
  //   1. Runs the base constructor (without setting up the problem).
  //   2. Sets up the problem.
  //      2.1. Converts the Eigen matrices to CSC matrices.
  //      2.2. Converts the vectors to dynamic arrays.
  //      2.3. Loads the problem formulation into the OSQP data object and sets up the workspace.
  //
  // Args:
  //   P: (n,n) matrix defining relations between parameters.
  //   A: (m,n) matrix defining parameter constraints relative to the lower and upper bound.
  //   q: (n) vector defining the linear cost of the problem.
  //   l: (m) vector defining the lower bound problem constraint.
  //   u: (m) vector defining the upper bound problem constraint.
  //   eps_abs: Absolute convergence tolerance.
  OSQPInterface(const Eigen::MatrixXd &P, const Eigen::MatrixXd &A, const std::vector<double> &q,
                const std::vector<double> &l, const std::vector<double> &u, const c_float eps_abs);

  // For freeing dynamic memory used by OSQP's data object.
  ~OSQPInterface();

  /****************
   * OPTIMIZATION
   ****************/
  // Solves the stored convec quadratic program (QP) problem using the OSQP solver.
  //
  // The function returns a tuple containing the solution as two float vectors.
  // The first element of the tuple contains the 'primal' solution. The second element contains the 'lagrange multiplier'
  // solution. The third element contains an integer with solver polish status information.
  //
  // How to use:
  //   1. Generate the Eigen matrices P, A and vectors q, l, u according to the problem.
  //   2. Initialize the interface and set up the problem.
  //        osqp_interface = OSQPInterface(P, A, q, l, u, 1e-6);
  //   3. Call the optimization function.
  //        std::tuple<std::vector<double>, std::vector<double>> result;
  //        result = osqp_interface.optimize();
  //   4. Access the optimized parameters.
  //        std::vector<float> param = std::get<0>(result);
  //        double x_0 = param[0];
  //        double x_1 = param[1];
  std::tuple<std::vector<double>, std::vector<double>, int> optimize();

  // Solves convex quadratic programs (QPs) using the OSQP solver.
  //
  // The function returns a tuple containing the solution as two float vectors.
  // The first element of the tuple contains the 'primal' solution. The second element contains the 'lagrange multiplier'
  // solution. The third element contains an integer with solver polish status information.
  //
  // How to use:
  //   1. Generate the Eigen matrices P, A and vectors q, l, u according to the problem.
  //   2. Initialize the interface.
  //        osqp_interface = OSQPInterface(1e-6);
  //   3. Call the optimization function with the problem formulation.
  //        std::tuple<std::vector<double>, std::vector<double>> result;
  //        result = osqp_interface.optimize(P, A, q, l, u, 1e-6);
  //   4. Access the optimized parameters.
  //        std::vector<float> param = std::get<0>(result);
  //        double x_0 = param[0];
  //        double x_1 = param[1];
  std::tuple<std::vector<double>, std::vector<double>, int> optimize(const Eigen::MatrixXd P,
                                                                const Eigen::MatrixXd A,
                                                                const std::vector<double> q,
                                                                const std::vector<double> l,
                                                                const std::vector<double> u);

  /**************************
   * DATA-RELATED FUNCTIONS
   **************************/

  // Converts the input data and sets up the workspace object.
  //
  // Args:
  //   P: (n,n) matrix defining relations between parameters.
  //   A: (m,n) matrix defining parameter constraints relative to the lower and upper bound.
  //   q: (n) vector defining the linear cost of the problem.
  //   l: (m) vector defining the lower bound problem constraint.
  //   u: (m) vector defining the upper bound problem constraint.
  c_int initializeProblem(const Eigen::MatrixXd &P,
                          const Eigen::MatrixXd &A,
                          const std::vector<double> &q,
                          const std::vector<double> &l,
                          const std::vector<double> &u);

  // Updates problem parameters while keeping solution in memory.
  //
  // Args:
  //   P_new: (n,n) matrix defining relations between parameters.
  //   A_new: (m,n) matrix defining parameter constraints relative to the lower and upper bound.
  //   q_new: (n) vector defining the linear cost of the problem.
  //   l_new: (m) vector defining the lower bound problem constraint.
  //   u_new: (m) vector defining the upper bound problem constraint.
  c_int updateP(const Eigen::MatrixXd &P_new);
  c_int updateA(const Eigen::MatrixXd &A_new);
  c_int updateQ(const std::vector<double> &q_new);
  c_int updateL(const std::vector<double> &l_new);
  c_int updateU(const std::vector<double> &u_new);
};

} // namespace osqp

#endif // OSQP_INTERFACE_H
