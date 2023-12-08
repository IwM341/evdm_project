cmake_minimum_required(VERSION 3.4)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")

option(VERBOSE "" ON)
if(VERBOSE)
   add_compile_definitions(_VERBOSE)
endif()

project(evdm_python_mudule LANGUAGES CXX)

find_package(pybind11 CONFIG REQUIRED)
message("byind found ${pybind11_FOUND}")

set(P_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/evdm_python)
set(SOURCES 
    ${P_SOURCES_DIR}/pyevdm.cpp
    ${P_SOURCES_DIR}/body_model.cpp
    ${P_SOURCES_DIR}/el_grid.cpp
    ${P_SOURCES_DIR}/core_python_defs.cpp
)

find_package(Eigen3 REQUIRED)

add_library(pyevdm MODULE ${SOURCES})
set_property(TARGET pyevdm PROPERTY CXX_STANDARD 17)

target_link_libraries(pyevdm pybind11::module pybind11::lto pybind11::windows_extras)
target_link_libraries(pyevdm Eigen3::Eigen)

target_include_directories(pyevdm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set(GROB_INCLIDE_DIR "D:/Tmp/GridObject_1/include")
target_include_directories(pyevdm PRIVATE ${GROB_INCLIDE_DIR})
pybind11_extension(pyevdm)

set(EIGEN_USE_MKL OFF)
if(EIGEN_USE_MKL)
    find_package()
    add_compile_definitions(EIGEN_USE_MKL_ALL)
endif()
