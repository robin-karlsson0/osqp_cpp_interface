# osqp_cpp_interface
Native C++ interface for the OSQP solver

## Setup the repository

1. Clone from Gitlab

2. Add the following folder structure
- osqp_cpp_interface/lib/osqp/lib
- osqp_cpp_interface/lib/osqp/include
- osqp_cpp_interface/build

## Install OSQP by building from souce

https://osqp.org/docs/get_started/sources.html

Following these instructions will provide the following library file and include files.
- osqp/build/out/libosqp.a
- osqp/include

The library file 'libosqp.a' together with the 'include' files will be used to run the solver.

1. Copy 'libosqp.a' from 'osqp' to the 'osqp_test' repository

osqp_cpp_interface/lib/osqp/lib/

2. Copy the 'include' files from 'osqp' to the 'osqp_test' repository

osqp_cpp_interface/lib/osqp/include/

## Build the test programs

Starting in the root:
1. cd build
2. cmake ../
3. make all

## Run the test programs

Starting in build/ directory:
- ./osqp_cpp_interface
- ./runTests

Successfully running 'osqp_cpp_interface' will provide the following output

-----------------------------------------------------------------
           OSQP v0.6.0  -  Operator Splitting QP Solver
              (c) Bartolomeo Stellato,  Goran Banjac
        University of Oxford  -  Stanford University 2019
-----------------------------------------------------------------
problem:  variables n = 2, constraints m = 3
          nnz(P) + nnz(A) = 9
settings: linear system solver = qdldl,
          eps_abs = 1.0e-03, eps_rel = 1.0e-03,
          eps_prim_inf = 1.0e-04, eps_dual_inf = 1.0e-04,
          rho = 1.00e-01 (adaptive),
          sigma = 1.00e-06, alpha = 1.00, max_iter = 4000
          check_termination: on (interval 25),
          scaling: on, scaled_termination: off
          warm start: on, polish: off, time_limit: off

iter   objective    pri res    dua res    rho        time
   1  -9.3904e+00   3.61e-01   3.18e+00   1.00e-01   6.83e-05s
  25  -8.2241e+00   3.08e-03   1.12e-03   1.00e-01   1.12e-04s

status:               solved
number of iterations: 25
optimal objective:    -8.2241
run time:             1.30e-04s
optimal rho estimate: 2.03e-01

0.668395
1.33266
