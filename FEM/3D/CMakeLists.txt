cmake_minimum_required(VERSION 3.5)
project(pybind_test VERSION 0.1.0)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(pybind11)

add_library(other MyClass_hexahedral_edge_element/node.cpp MyClass_hexahedral_edge_element/element.cpp)
set_target_properties(other PROPERTIES LINKER_LANGUAGE CXX)

pybind11_add_module(FEM_3D FEM_3D.cpp)

target_link_libraries(FEM_3D PRIVATE other)

# otherライブラリのincludeパスを設定
target_include_directories(other PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/MyClass_hexahedral_edge_element
)

# FEM_3Dモジュールのincludeパスも設定
target_include_directories(FEM_3D PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/MyClass_hexahedral_edge_element
)