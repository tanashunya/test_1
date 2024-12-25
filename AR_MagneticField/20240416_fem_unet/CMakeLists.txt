cmake_minimum_required(VERSION 3.5)
project(pybind_test VERSION 0.1.0)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(pybind11)

add_library(other GenMatrix.cpp MyFunction.cpp)

pybind11_add_module(Carray MyModule.cpp)
pybind11_add_module(Cfunctions MyModule2.cpp)

target_link_libraries(Carray PRIVATE other)
target_link_libraries(Cfunctions PRIVATE other)