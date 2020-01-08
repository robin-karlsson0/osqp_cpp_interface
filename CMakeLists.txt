cmake_minimum_required(VERSION 2.8.3)
project(osqp_cpp_interface)

add_compile_options(-std=c++11)

enable_testing()

include_directories(
    include
    lib/osqp/include
)

add_executable(osqp_cpp_interface
    src/main.cpp
    src/osqp_interface.cpp
    src/csc_matrix_conv.cpp
)

target_link_libraries(osqp_cpp_interface
    ${CMAKE_SOURCE_DIR}/lib/osqp/lib/libosqp.a
    ${CMAKE_DL_LIBS}
)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(runTests
  test/src/test_osqp_interface.cpp
  test/src/test_csc_matrix_conv.cpp
  src/osqp_interface.cpp
  src/csc_matrix_conv.cpp
)

target_link_libraries(runTests
    ${CMAKE_SOURCE_DIR}/lib/osqp/lib/libosqp.a
    ${CMAKE_DL_LIBS}
)

target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)
