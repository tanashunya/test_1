cmake_minimum_required(VERSION 3.5)
project(pybind_test VERSION 0.1.0)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(pybind11)

add_library(other GenMatrix.cpp MyFunction.cpp)

pybind11_add_module(MOTOR_FEM Motor_FEM.cpp)

target_link_libraries(MOTOR_FEM PRIVATE other)