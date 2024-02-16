cmake_minimum_required(VERSION 3.4)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")


option(VERBOSE "" ON)
if(VERBOSE)
   add_compile_definitions(_VERBOSE)
endif()

#types, available in structs. 
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/type_properties.cmake)
#optional



project(evdm_python_mudule LANGUAGES CXX VERSION 1.2.3)

find_package(pybind11 CONFIG REQUIRED)
message("pybind found ${pybind11_FOUND}")

set(P_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/evdm_python)
set(FORM_FACTOR_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/NuclearFormFactors)
set(SOURCES 
    ${P_SOURCES_DIR}/pyevdm.cpp
    ${P_SOURCES_DIR}/body_model.cpp
    ${P_SOURCES_DIR}/el_grid.cpp
    ${P_SOURCES_DIR}/core_python_defs.cpp
    ${P_SOURCES_DIR}/distribution.cpp
    ${P_SOURCES_DIR}/capture.cpp
    ${P_SOURCES_DIR}/matrix.cpp
    ${P_SOURCES_DIR}/pyscatter.cpp
)

include(cmake/CPM.cmake)

#CPMAddPackage("gh:Iwm341/GridObjects")


find_package(Eigen3 REQUIRED)

add_library(pyevdm MODULE ${SOURCES})
set_property(TARGET pyevdm PROPERTY CXX_STANDARD 20)

target_link_libraries(pyevdm pybind11::module pybind11::lto pybind11::windows_extras)
target_link_libraries(pyevdm Eigen3::Eigen)

target_include_directories(pyevdm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

## grob lib, will be downloaded with cpm_cmake
set(GROB_INCLIDE_DIR "D:/Tmp/GridObject_1/include")
target_include_directories(pyevdm PRIVATE ${GROB_INCLIDE_DIR})
##

pybind11_extension(pyevdm)

set(EIGEN_USE_MKL OFF)
if(EIGEN_USE_MKL)
    find_package()
    add_compile_definitions(EIGEN_USE_MKL_ALL)
endif()



#copy python module to output file
set(PYMODULE_DIR ${CMAKE_BINARY_DIR}/evdm)
set(PYMODULE_SRC_DIR ${CMAKE_SOURCE_DIR}/evdm_python_package_files)
file(COPY ${PYMODULE_SRC_DIR}/evdm DESTINATION ${PYMODULE_DIR})
message(PROJECT_VERSION = ${PROJECT_VERSION})
configure_file(${PYMODULE_SRC_DIR}/pyproject.toml.in ${PYMODULE_DIR}/pyproject.toml)
#
