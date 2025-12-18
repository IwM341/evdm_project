cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)


project(evdm_python_mudule LANGUAGES CXX VERSION 2.0.0)

#include("build/build/generators/conan_toolchain.cmake")

option(GEN_R_NO_POW "in capture r generator, use r = g() instead g()^p" OFF)
if(${GEN_R_NO_POW})
   add_compile_definitions(NO_POW)
endif()

option(VERBOSE "Debug some functions" OFF)
if(${VERBOSE})
   add_compile_definitions(_VERBOSE)
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
option(OPENMP "" ON)
option(CATCHNAN "" ON)
else()
option(OPENMP "" OFF)
option(CATCHNAN "" OFF)
endif()

if(${CATCHNAN})
   add_compile_definitions(CATCHNAN)
endif()

option(BUILD_TESTS "BUILD_TESTS" OFF)

message("USEOPENMP: ${OPENMP}" )
message("CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE})
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug flags = " ${CMAKE_CXX_FLAGS})
endif()
if(CMAKE_BUILD_TYPE MATCHES Release)
    message("Release flags = " ${CMAKE_CXX_FLAGS})
    message("Release flags = " ${CMAKE_CXX_FLAGS_RELEASE})
endif()
if(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    message("RelWithDebInfo flags = " ${CMAKE_CXX_FLAGS})
endif()

#types, available in structs. 
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/type_properties.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/scatter_config.cmake)
#optional







if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  message(Clang)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  message(GCC)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(Intel)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(MSVC)
else()
    message(ERROR)
endif()

if( 1 OR CMAKE_BUILD_TYPE MATCHES Debug OR
       CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    
    if (MSVC)
        add_compile_options(/bigobj)
        message("/bigobj option set")
    else ()
        #add_compile_options(-Wa,-mbig-obj)
        message("-Wa,-mbig-obj option set")
    endif()
endif()

set(P_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/evdm_python)
set(FORM_FACTOR_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/NuclearFormFactors)

file(GLOB_RECURSE SOURCES "${P_SOURCES_DIR}/*.cpp")



include(cmake/CPM.cmake)

find_package(Eigen3)
if(NOT Eigen3_FOUND)
    message("Eigen3 not found, download from git")
    CPMAddPackage(NAME Eigen GIT_REPOSITORY  "https://gitlab.com/libeigen/eigen" GIT_TAG 3.4.0)
endif()
find_package(pybind11 CONFIG)
if(NOT pybind11_FOUND)
    message("pybind11 not found, download from git")
    CPMAddPackage("gh:pybind/pybind11@2.13.6")
endif()


find_package(Grob)
if(NOT Grob_FOUND )
    message("Grob not found, download from git")
    CPMAddPackage("gh:IwM341/GridObject@1.1.4")
endif()




add_library(pyevdm MODULE ${SOURCES} )
set_property(TARGET pyevdm PROPERTY CXX_STANDARD 20)

target_link_libraries(pyevdm pybind11::module)
target_link_libraries(pyevdm Eigen3::Eigen)



message("boost dir ${Boost_INCLUDE_DIR}")

if(OPENMP)
    message("use openmp")
    find_package(OpenMP REQUIRED)
    target_link_libraries(pyevdm OpenMP::OpenMP_CXX)
endif()



target_include_directories(pyevdm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

## grob lib, will be downloaded with cpm_cmake
#target_include_directories(pyevdm PRIVATE ${GROB_INCLUDE_DIR})
#target_include_directories(pyevdm PRIVATE "D:/Tmp/GridObject_1/include")
target_link_libraries(pyevdm Grob::grob)
##

pybind11_extension(pyevdm)





##copy python module to output file
set(PYMODULE_DIR ${CMAKE_BINARY_DIR}/evdm)
set(PYMODULE_SRC_DIR ${CMAKE_SOURCE_DIR}/evdm_python_package_files)
file(COPY ${PYMODULE_SRC_DIR}/evdm DESTINATION ${PYMODULE_DIR})
#file(COPY ${PYMODULE_SRC_DIR}/evdm-stubs DESTINATION ${PYMODULE_DIR})


##after development delete
file(REMOVE ${PYMODULE_DIR}/evdm/_cpp_lib/pyevdm.py)
##

message(PROJECT_VERSION = ${PROJECT_VERSION})
configure_file(${PYMODULE_SRC_DIR}/pyproject.toml.in ${PYMODULE_DIR}/pyproject.toml)

set_target_properties(pyevdm PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PYMODULE_DIR}/evdm/_cpp_lib
    LIBRARY_OUTPUT_DIRECTORY ${PYMODULE_DIR}/evdm/_cpp_lib
)



set (TEST_DIR ${CMAKE_SOURCE_DIR}/tests)
message("BUILD_TESTS = ${BUILD_TESTS}")

function (prepare_test TESTNAME TESTPATH)
    add_executable(${TESTNAME} ${TESTPATH})
    set_property(TARGET ${TESTNAME} PROPERTY CXX_STANDARD 20)
    message("test = " ${TESTNAME})
    message("file = " ${TESTPATH})
    target_link_libraries(${TESTNAME} Grob::grob)
    target_link_libraries(${TESTNAME} Eigen3::Eigen)
    if(OPENMP)
        target_link_libraries(${TESTNAME} OpenMP::OpenMP_CXX)
    endif()
    target_include_directories(${TESTNAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${TESTNAME} PRIVATE ${GROB_INCLUDE_DIR})
endfunction()

if(${BUILD_TESTS})
    file(GLOB files "${TEST_DIR}/*.cpp")
    foreach(filename ${files})
        get_filename_component(TEST_NAME ${filename} NAME_WE)
        set(test_no_pow  "${TEST_NAME}_nopow") 

        prepare_test(${test_no_pow} ${filename})
        target_compile_definitions(${test_no_pow} PRIVATE NO_POW)

        prepare_test(${TEST_NAME} ${filename})
    
    endforeach() 
endif()
