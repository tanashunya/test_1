cmake_minimum_required(VERSION 3.5)
project(pybind_test VERSION 0.1.0)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(pybind11)

add_library(other node.cpp element.cpp edge.cpp)
set_target_properties(other PROPERTIES LINKER_LANGUAGE CXX)

pybind11_add_module(FEM_3D FEM_3D.cpp)

target_link_libraries(FEM_3D PRIVATE other)
