#include "osqp_interface.h"
#include <iostream>

namespace osqp
{

OSQPInterface::OSQPInterface(const c_float eps_abs)
{
  /************************
   * INITIALIZE WORKSPACE
   ************************/
  settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  data = (OSQPData *)c_malloc(sizeof(OSQPData));

  /*******************
   * SOLVER SETTINGS
   *******************/
  if (settings)
  {
    osqp_set_default_settings(settings);
    settings->alpha = 1.0;
    settings->eps_abs = eps_abs;
  }
  // Set flag for successful initialization
  exitflag = 0;
}

OSQPInterface::OSQPInterface(const Eigen::MatrixXd &P, const Eigen::MatrixXd &A, const std::vector<double> &q,
                             const std::vector<double> &l, const std::vector<double> &u, const c_float eps_abs) : OSQPInterface(eps_abs)
{
  initializeProblem(P, A, q, l, u);
}

c_int OSQPInterface::initializeProblem(const Eigen::MatrixXd &P,
                                       const Eigen::MatrixXd &A,
                                       const std::vector<double> &q,
                                       const std::vector<double> &l,
                                       const std::vector<double> &u)
{
  /*******************
   * SET UP MATRICES
   *******************/
  // Transform 'P' into an 'upper trapesoidal matrix'
  Eigen::MatrixXd P_trap = P.triangularView<Eigen::Upper>();
  // CSC matrices
  CSC_Matrix P_csc = convEigenMatrixToCSCMatrix(P_trap);
  CSC_Matrix A_csc = convEigenMatrixToCSCMatrix(A);
  // Dynamic float arrays
  std::vector<double> q_tmp(q.begin(), q.end());
  std::vector<double> l_tmp(l.begin(), l.end());
  std::vector<double> u_tmp(u.begin(), u.end());
  double *q_dyn = q_tmp.data();
  double *l_dyn = l_tmp.data();
  double *u_dyn = u_tmp.data();
  //int *element_n;
  //CSC_Matrix P_csc = transformP(P, element_n);
  //CSC_Matrix A_csc = transformA(A);
  //double* q_dyn = transformQ(q);
  //double* l_dyn = transformL(l);
  //double* u_dyn = transformU(u);

  /**********************
   * OBJECTIVE FUNCTION
   **********************/
  // Number of constraints
  c_int constr_m = A.rows();
  // Number of parameters
  param_n = P.rows();
  // Number of elements
  c_int P_elem_N = P_trap.nonZeros();
  c_int A_elem_N = A.nonZeros();

  /*****************
   * POLULATE DATA
   *****************/
  data->m = constr_m;
  data->n = param_n;
  data->P = csc_matrix(data->n, data->n, P_elem_N, P_csc.elem_val.data(), P_csc.row_idx.data(), P_csc.col_idx.data());
  data->q = q_dyn;
  data->A = csc_matrix(data->m, data->n, A_elem_N, A_csc.elem_val.data(), A_csc.row_idx.data(), A_csc.col_idx.data());
  data->l = l_dyn;
  data->u = u_dyn;

  // For deconstructor
  problem_in_memory = true;

  // Setup workspace
  exitflag = osqp_setup(&work, data, settings);

  return exitflag;
}

OSQPInterface::~OSQPInterface()
{
  // Cleanup dynamic OSQP memory
  if (data)
  {
    if (problem_in_memory)
    {
      c_free(data->A);
      c_free(data->P);
    }
    c_free(data);
  }
  if (settings)
    c_free(settings);
}

c_int OSQPInterface::updateP(const Eigen::MatrixXd &P_new)
{
  // Transform 'P' into an 'upper trapesoidal matrix'
  Eigen::MatrixXd P_new_trap = P_new.triangularView<Eigen::Upper>();
  // CSC matrices
  CSC_Matrix P_csc = convEigenMatrixToCSCMatrix(P_new_trap);
  // Number of elements
  c_int P_elem_N = P_new_trap.nonZeros();
  osqp_update_P(work, P_csc.elem_val.data(), OSQP_NULL, P_elem_N);
}

c_int OSQPInterface::updateA(const Eigen::MatrixXd &A_new)
{
  // CSC matrices
  CSC_Matrix A_csc = convEigenMatrixToCSCMatrix(A_new);
  // Number of elements
  c_int A_elem_N = A_new.nonZeros();
  osqp_update_A(work, A_csc.elem_val.data(), OSQP_NULL, A_elem_N);
}

c_int OSQPInterface::updateQ(const std::vector<double> &q_new)
{
  std::vector<double> q_tmp(q_new.begin(), q_new.end());
  double *q_dyn = q_tmp.data();
  osqp_update_lin_cost(work, q_dyn);
}

c_int OSQPInterface::updateL(const std::vector<double> &l_new)
{
  std::vector<double> l_tmp(l_new.begin(), l_new.end());
  double *l_dyn = l_tmp.data();
  osqp_update_lower_bound(work, l_dyn);
}

c_int OSQPInterface::updateU(const std::vector<double> &u_new)
{
  std::vector<double> u_tmp(u_new.begin(), u_new.end());
  double *u_dyn = u_tmp.data();
  osqp_update_upper_bound(work, u_dyn);
}

std::tuple<std::vector<double>, std::vector<double>> OSQPInterface::solve()
{
  // Solve Problem
  osqp_solve(work);

  /********************
   * EXTRACT SOLUTION
   ********************/
  double *sol_x = work->solution->x;
  double *sol_y = work->solution->y;
  std::vector<double> sol_primal(sol_x, sol_x + static_cast<int>(param_n));
  std::vector<double> sol_lagrange_multiplier(sol_y, sol_y + static_cast<int>(param_n));
  // Result tuple
  std::tuple<std::vector<double>, std::vector<double>> result = std::make_tuple(sol_primal, sol_lagrange_multiplier);

  return result;
}

std::tuple<std::vector<double>, std::vector<double>> OSQPInterface::optimize()
{
  // Run the solver on the stored problem representation.
  std::tuple<std::vector<double>, std::vector<double>> result = solve();
  return result;
}

std::tuple<std::vector<double>, std::vector<double>> OSQPInterface::optimize(const Eigen::MatrixXd P,
                                                                             const Eigen::MatrixXd A,
                                                                             const std::vector<double> q,
                                                                             const std::vector<double> l,
                                                                             const std::vector<double> u)
{
  initializeProblem(P, A, q, l, u);

  // Run the solver on the stored problem representation.
  std::tuple<std::vector<double>, std::vector<double>> result = solve();
  return result;
}

inline bool OSQPInterface::isEqual(double x, double y)
{
  const double epsilon = 1e-6;
  return std::abs(x - y) <= epsilon * std::abs(x);
}

CSC_Matrix OSQPInterface::convEigenMatrixToCSCMatrix(const Eigen::MatrixXd A)
{
  // Input dense matrix dimensions
  int A_rows = A.rows();
  int A_cols = A.cols();

  /************************************************************************
   * Generate 'sparse matrix B' from nonzero elements in 'dense matrix A'
   ************************************************************************/

  // Generate list of nonzero elements
  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(A.size());
  for (int i = 0; i < A_rows; i++)
  {
    for (int j = 0; j < A_cols; j++)
    {
      double A_val = A(i, j);
      if (!isEqual(A_val, 0.0))
      {
        triplet_list.push_back(Eigen::Triplet<double>(i, j, A_val));
      }
    }
  }
  // Generate 'sparse matrix B' and fill with nonzero elements in list
  Eigen::SparseMatrix<double> B(A_rows, A_cols);
  B.setFromTriplets(triplet_list.begin(), triplet_list.end());

  /************************************************************************
   * Generate 'Compressed Sparse Column (CSC) Matrix A_csc' struct object
   ************************************************************************/

  // Generate CSC matrix B
  int B_nonzero_N = B.nonZeros();
  B.makeCompressed();

  // Extract pointer arrays
  double *val_ptr = B.valuePtr();
  int *inn_ptr = B.innerIndexPtr();
  int *out_ptr = B.outerIndexPtr();

  // Copy values of pointer arrays into vectors in CSC struct
  // Array lengths:
  //     elem_val : nonzero element count
  //     row_idx  : nonzero element count
  //     col_idx  : input matrix column count + 1
  CSC_Matrix A_csc;
  A_csc.elem_val.assign(val_ptr, val_ptr + B_nonzero_N);
  A_csc.col_idx.assign(out_ptr, out_ptr + A_cols + 1);
  A_csc.row_idx.assign(inn_ptr, inn_ptr + B_nonzero_N);

  return A_csc;
}

c_int OSQPInterface::getExitFlag(void)
{
  return exitflag;
}

/*
CSC_Matrix OSQPInterface::transformP(const Eigen::MatrixXd &P, int *nonzeros)
{
  // Transform 'P' into an 'upper trapesoidal matrix'
  Eigen::MatrixXd P_trap = P.triangularView<Eigen::Upper>();
  *nonzeros = P_trap.nonZeros();
  // CSC matrices
  CSC_Matrix P_csc = convEigenMatrixToCSCMatrix(P_trap);
  return P_csc;
}

CSC_Matrix OSQPInterface::transformA(const Eigen::MatrixXd &A)
{
  CSC_Matrix A_csc = convEigenMatrixToCSCMatrix(A);
  return A_csc;
}

double* OSQPInterface::transformQ(const std::vector<double> &q)
{
  // Dynamic float arrays
  std::vector<double> q_tmp(q.begin(), q.end());
  double* q_dyn = q_tmp.data();
  return q_dyn;
}

double* OSQPInterface::transformL(const std::vector<double> &l)
{
  std::vector<double> l_tmp(l.begin(), l.end());
  double* l_dyn = l_tmp.data();
  return l_dyn;
}

double* OSQPInterface::transformU(const std::vector<double> &u)
{
  std::vector<double> u_tmp(u.begin(), u.end());
  double* u_dyn = u_tmp.data();
  return u_dyn;
}
*/

} // namespace osqp
