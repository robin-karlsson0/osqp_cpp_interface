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
  // Transform 'P' into a sparse matrix and extract data as dynamic arrays
  Eigen::SparseMatrix<double> P_sparse = P_trap.sparseView();
  double *P_val_ptr = P_sparse.valuePtr();
  int *P_row_ptr_tmp = P_sparse.innerIndexPtr();
  int *P_col_ptr_tmp = P_sparse.outerIndexPtr();
  // Convert dynamic 'int' arrays to 'c_int' arrays (OSQP input type)
  c_int P_elem_N = P_sparse.nonZeros();
  c_int *P_row_ptr = new c_int[P_elem_N];
  for (int i = 0; i < P_elem_N; i++)
  {
    P_row_ptr[i] = c_int(P_row_ptr_tmp[i]);
  }
  c_int *P_col_ptr = new c_int[P_elem_N];
  for (int i = 0; i < P_elem_N; i++)
  {
    P_col_ptr[i] = c_int(P_col_ptr_tmp[i]);
  }

  // Transform 'A' into a sparse matrix and extract data as dynamic arrays
  Eigen::SparseMatrix<double> A_sparse = A.sparseView();
  double *A_val_ptr = A_sparse.valuePtr();
  int *A_row_ptr_tmp = A_sparse.innerIndexPtr();
  int *A_col_ptr_tmp = A_sparse.outerIndexPtr();
  // Convert dynamic 'int' arrays to 'c_int' arrays (OSQP input type)
  c_int A_elem_N = A_sparse.nonZeros();
  c_int *A_row_ptr = new c_int[A_elem_N];
  for (int i = 0; i < A_elem_N; i++)
  {
    A_row_ptr[i] = c_int(A_row_ptr_tmp[i]);
  }
  c_int *A_col_ptr = new c_int[A_elem_N];
  for (int i = 0; i < A_elem_N; i++)
  {
    A_col_ptr[i] = c_int(A_col_ptr_tmp[i]);
  }
  // Dynamic float arrays
  std::vector<double> q_tmp(q.begin(), q.end());
  std::vector<double> l_tmp(l.begin(), l.end());
  std::vector<double> u_tmp(u.begin(), u.end());
  double *q_dyn = q_tmp.data();
  double *l_dyn = l_tmp.data();
  double *u_dyn = u_tmp.data();

  /**********************
   * OBJECTIVE FUNCTION
   **********************/
  // Number of constraints
  c_int constr_m = A.rows();
  // Number of parameters
  param_n = P.rows();

  /*****************
   * POLULATE DATA
   *****************/
  data->m = constr_m;
  data->n = param_n;
  data->P = csc_matrix(data->n, data->n, P_elem_N, P_val_ptr, P_row_ptr, P_col_ptr);
  data->q = q_dyn;
  data->A = csc_matrix(data->m, data->n, A_elem_N, A_val_ptr, A_row_ptr, A_col_ptr);
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
  Eigen::MatrixXd P_trap = P_new.triangularView<Eigen::Upper>();
  // Transform 'P' into a sparse matrix and extract data as dynamic arrays
  Eigen::SparseMatrix<double> P_sparse = P_trap.sparseView();
  double *P_val_ptr = P_sparse.valuePtr();
  // Convert dynamic 'int' arrays to 'c_int' arrays (OSQP input type)
  c_int P_elem_N = P_sparse.nonZeros();
  osqp_update_P(work, P_val_ptr, OSQP_NULL, P_elem_N);
}

c_int OSQPInterface::updateA(const Eigen::MatrixXd &A_new)
{
  // Transform 'A' into a sparse matrix and extract data as dynamic arrays
  Eigen::SparseMatrix<double> A_sparse = A_new.sparseView();
  double *A_val_ptr = A_sparse.valuePtr();
  // Convert dynamic 'int' arrays to 'c_int' arrays (OSQP input type)
  c_int A_elem_N = A_sparse.nonZeros();

  osqp_update_A(work, A_val_ptr, OSQP_NULL, A_elem_N);
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

std::tuple<std::vector<double>, std::vector<double>, int> OSQPInterface::solve()
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
  // Solver polish status
  int status_polish = work->info->status_polish;
  // Result tuple
  std::tuple<std::vector<double>, std::vector<double>, int> result = std::make_tuple(sol_primal, sol_lagrange_multiplier, status_polish);

  return result;
}

std::tuple<std::vector<double>, std::vector<double>, int> OSQPInterface::optimize()
{
  // Run the solver on the stored problem representation.
  std::tuple<std::vector<double>, std::vector<double>, int> result = solve();
  return result;
}

std::tuple<std::vector<double>, std::vector<double>, int> OSQPInterface::optimize(const Eigen::MatrixXd P,
                                                                                  const Eigen::MatrixXd A,
                                                                                  const std::vector<double> q,
                                                                                  const std::vector<double> l,
                                                                                  const std::vector<double> u)
{
  initializeProblem(P, A, q, l, u);

  // Run the solver on the stored problem representation.
  std::tuple<std::vector<double>, std::vector<double>, int> result = solve();
  return result;
}

inline bool OSQPInterface::isEqual(double x, double y)
{
  const double epsilon = 1e-6;
  return std::abs(x - y) <= epsilon * std::abs(x);
}

c_int OSQPInterface::getExitFlag(void)
{
  return exitflag;
}

} // namespace osqp
